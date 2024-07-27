from multiprocessing import Process, Queue,Lock
# import threading
import os

from PySide6.QtWidgets import QTextEdit,QVBoxLayout, QWidget,QListWidget,QLabel,QGridLayout,QHBoxLayout,QListWidgetItem
from PySide6.QtCore import QTimer, Signal, QObject,Qt
from PySide6.QtGui import QTextCursor,QColor

from pipeline import cal_library

class LogOutput(QTextEdit):
    def __init__(self, parent=None):
        super().__init__(parent)
        self.setReadOnly(True)
        self.setStyleSheet(
            """
            QTextEdit {
                background-color: transparent;
                color: #000000;
                font-family: Consolas, monospace;
                font-size: 16px;
                border: none;
            }
            """
        )

    def write(self, message):
        cursor = self.textCursor()
        cursor.movePosition(QTextCursor.End)
        cursor.insertText(message)
        self.ensureCursorVisible()
        self.verticalScrollBar().setValue(self.verticalScrollBar().maximum())

    def on_scoll(self):
        pass

    def flush(self):
        self.clear()


class Worker(QObject):
    success_signal = Signal()  # 参数设置成功
    fail_signal= Signal() # 参数设置失败
    error_signal=Signal() # 出现错误
    
    lib_signal=Signal() # lib 计算标识
    
    def __init__(
        self,
        id: int,
        function,
        func_params: dict,
        mzxml_file_id: int,
        log_output: LogOutput,
        parent=None
    ):
        """
        TODO: 自行填补所需的参数, 我这边列举了一些参数，你可以根据自己的设计方案去增删改查参数
        -   function: 执行的目标函数
        -   func_params: 函数参数 (因为除了文件名以外其他参数是一样的)
        -   mzxml_file_id: 需要分析的文件下标
        -   LogOutput(QTextEdit): 用于显示展示日志的页面 
        """
        super().__init__(parent)
        self.id = id
        self.func = function
        self.func_params = func_params
        self.log_output = log_output
        self.mzxml_file_id = mzxml_file_id

        # 定时器触发执行操作 (将执行到哪一步的信息发给 view_interface)
        self.text_queue = Queue()
        self.timer = QTimer()
        self.process = None
        self._connect_signal2solt()

    def _connect_signal2solt(self):
        """
        TODO: 将 process 结束的信号连接到 end_signal 上, 最好是用根据 log_output 输出 'over' 作为标志
        """
        self.timer.timeout.connect(self._check_queue)
        # self.end_signal.connect()
        # self.process.finished.connect(self.end_signal)

    def _check_queue(self):
        if not self.text_queue.empty():
            result = self.text_queue.get()
            self.log_output.write(result)
            if result == 'over success':
                self.success_signal.emit()
            elif result == 'over failing':
                self.fail_signal.emit()
            elif result == 'over error':
                self.error_signal.emit()
            elif result == 'cal library over':
                self.lib_signal.emit()

    def run(self):
        """
        TODO: 执行 function 函数, 启动多进程/启动计时器
        """
        kwargs = {**self.func_params, 'text_queue': self.text_queue}
        self.process = Process(target=self.func,kwargs=kwargs)
        self.timer.start(100)
        self.process.start()

    def destroy_resource(self):
        """
        TODO: 释放进程资源，销毁 worker
        """
        if self.process is not None and self.process.is_alive():
            self.process.terminate()
            self.process.join()
        if self.timer:
            self.timer.stop()

    def reset_run(self, func, func_params, mzxml_file_id: int):
        del self.text_queue
        self.text_queue = Queue()
        self.log_output.flush()
        self.func = func
        self.func_params = func_params
        self.mzxml_file_id = mzxml_file_id

    def reset_worker(self, func, func_params, mzxml_file_id):
        self.destroy_resource()
        self.reset_run(func, func_params, mzxml_file_id)
        self.run()
        
            
    def callib_first(self):
        # cal_library
        kwargs = {**self.func_params, 'text_queue': self.text_queue}
        self.process = Process(target=cal_library,kwargs=kwargs)
        self.timer.start(100)
        self.process.start()
        
        


class ViewInterface(QWidget):
    over_signal = Signal()

    def __init__(
        self,
        object_name: str,
        parent=None
    ):
        """
        TODO: 增加 worker_num 参数，用于控制并行化的进程数量

        可以添加一个 worker 列表表示目前分析的 worker 数量 (查看目前正在分析的文件)
        """
        super().__init__(parent)
        self.worker_num = 1 
        self.setObjectName(object_name)
        self.workers: dict[int,Worker] = {}
        self.lock=Lock()
        self.log_output=LogOutput(self)

       
        self.worker_logs = {}  # 存储每个 worker 的日志框

        layout = QVBoxLayout(self)

        self.setObjectName(object_name)

        # 创建文件状态列表
        files_layout=QHBoxLayout()
        self.unprocessed_list = QListWidget(self)
        self.processing_list = QListWidget(self)
        self.processed_list = QListWidget(self)
        files_layout.addWidget(QLabel("未处理的文件"))
        files_layout.addWidget(self.unprocessed_list)
        files_layout.addWidget(QLabel("正在处理的文件"))
        files_layout.addWidget(self.processing_list)
        files_layout.addWidget(QLabel("已处理完成的文件"))
        files_layout.addWidget(self.processed_list)
        layout.addLayout(files_layout)
        
        
        # 添加 Worker 日志框
        self.worker_log_layout = QGridLayout()
        layout.addLayout(self.worker_log_layout)

        
    def create_worker_log(self, worker_id):
        log_output = LogOutput(self)

        # 计算行列索引
        row = len(self.worker_logs) // 2  # 每行2个
        col = len(self.worker_logs) % 2   # 当前列数

        self.worker_log_layout.addWidget(QLabel(f"Worker {worker_id}"), row, col * 2)  # 标签
        self.worker_log_layout.addWidget(log_output, row, col * 2 + 1)  # 日志框
        self.worker_logs[worker_id] = log_output  # 保存日志框的引用
    
    def clear_worker_log(self):
        for i in reversed(range(self.worker_log_layout.count())):
            item=self.worker_log_layout.itemAt(i)
            wigdet=item.widget()
            self.worker_log_layout.removeWidget(wigdet)
            wigdet.deleteLater()
        for worker_log in self.worker_logs.values():
            worker_log.clear()
            del worker_log
        self.worker_logs.clear()

        
    def worker_start(self):
        for id in range(self.worker_num):
            with self.lock :
                self.processing_files[id]=self.cur_index
                self.processing_list.clear()
                for key,value in self.processing_files.items():
                    filename='Finished!'
                    if value!=-1:
                        filename=self.filenames[value]
                    self.processing_list.addItem(f"[Worker {key}]: {filename}")
                    
                unprocessed_name=self.filenames[self.cur_index]
                row=self.unprocessed_list.row(self.unprocessed_list.findItems(unprocessed_name,Qt.MatchFlag.MatchExactly)[0])
                self.unprocessed_list.takeItem(row)
                    
                worker=self.workers[id]
                func_params=self.other_params.copy()
                self.extend_func_params(
                    func_params, self.mzxml_files[self.cur_index])
                worker.reset_worker(self.function,func_params,self.cur_index)
                self.cur_index += 1
        
        

    def _init_workers(self, function, worker_num, mzxml_files, other_params):
        """
        初始化 worker 

        Parameters
        ---
        -   function: 并行化的函数
        -   mzxml_files: mzxml 需要分析的文件列表
        -   kwargs:   函数对应的参数 (mzxml_files 这里被单独拎出来，因为对它进行调度分配)

        TODO: 

        根据 worker_num 的值，初始先将前 worker_num 个 mzxml_file 分配到 worker_num 个 worker
        """
        self.function = function
        self.mzxml_files = mzxml_files
        self.filenames = [os.path.basename(file) for file in mzxml_files]
        self.worker_num = worker_num
        self.other_params = other_params
        self.processed_files:dict[int,list] = {}  # 下标
        self.processing_files:dict[int,int] ={} # 下标
        self.cur_index = 0

        self.unprocessed_list.clear()
        self.processing_list.clear()
        self.processed_list.clear()
        self.unprocessed_list.addItems(self.filenames)
        
        self.clear_worker_log()
        
        for id in range(self.worker_num):
            self.processing_files[id]=-1
            self.processed_files[id]=[]
        for id in range(self.worker_num):
            self.create_worker_log(id)
        # self.callib_first()
        for id in range(self.worker_num):
            self.set_worker(id, self.function, self.other_params.copy())
            
        libworker=self.workers[0]
        libworker.lib_signal.connect(lambda: self.worker_start())
        libworker.callib_first()

        

    def extend_func_params(self,func_params: dict, mzxml_file):
        func_params.update({
            'mzxml_files': [mzxml_file]
        })

    def set_worker(self, id: int, function, func_params: dict):
        """
        TODO: 
        把这个函数功能更改为有一个 worker 结束之后，分配下一个 mzxml_file 重新分配给它

        这里可以选择重启这个进程 (需要自己构造函数实现这个功能)，或者把 worker 资源释放，重新实例化 worker 类再启动

        具体采用什么策略你可以自己设计，我这边只是给一些建议
        """
        self.extend_func_params(
            func_params, self.mzxml_files[self.cur_index])
        
        worker = Worker(
            id, function, func_params,
            self.cur_index, self.worker_logs[id],parent=None
        )
        worker.success_signal.connect(lambda: self.switch_worker(worker.id, function, self.other_params.copy(),'Success'))
        worker.fail_signal.connect(lambda: self.switch_worker(worker.id, function,self.other_params.copy(),'Fail'))
        worker.error_signal.connect(lambda: self.switch_worker(worker.id, function,self.other_params.copy(),'Error'))
        self.workers[id] = worker
        
        
       

    def switch_worker(self, id: int, func, func_params: dict,type:str):
        """
        Paramerters:
        ---
        id: 介于 0 ~ worker_num -1 之间, 用于获取 worker 的位置
        TODO: 添加锁, 建议直接用粗粒度锁，直接将整个操作锁起来
        """
        with self.lock:
            self.processed_files[id].append(self.processing_files[id])
            if type=='Success':
                message='[Success]: '+self.filenames[self.processing_files[id]]
                item=QListWidgetItem(message)
                item.setForeground(QColor("green"))
                self.processed_list.addItem(item)
            elif type=='Fail':
                message='[Fail]: '+self.filenames[self.processing_files[id]]
                item=QListWidgetItem(message)
                item.setForeground(QColor("blue"))
                self.processed_list.addItem(item)
            elif type=='Error':
                message='[Error]: '+self.filenames[self.processing_files[id]]
                item=QListWidgetItem(message)
                item.setForeground(QColor("red"))
                self.processed_list.addItem(item)
                
                
            self.processing_files[id]=-1
            
            if self.cur_index>=len(self.mzxml_files):
                self.processing_list.clear()
                for key,value in self.processing_files.items():
                    filename='Finished!'
                    if value!=-1:
                        filename=self.filenames[value]
                    self.processing_list.addItem(f"[Worker {key}]: {filename}")
                
                self.worker_logs[id].write(f'\nWorker {id} has finished all tasks.')
                fin_num=0
                for i in range(self.worker_num):
                    if self.processing_files[i]==-1:
                        fin_num+=1
                if fin_num==self.worker_num:
                    for i in range(self.worker_num):
                        self.worker_logs[i].flush()
                        self.worker_logs[i].write(f'\nAll processed has been terminated.')
                    self.task_over()
                    # self.over_signal.emit()
            else:
                self.extend_func_params(func_params, self.mzxml_files[self.cur_index])
                worker = self.workers[id]
                
                self.processing_files[id]=self.cur_index
                self.processing_list.clear()
                for key,value in self.processing_files.items():
                    filename='Finished!'
                    if value!=-1:
                        filename=self.filenames[value]
                    self.processing_list.addItem(f"[Worker {key}]: {filename}")
                
                unprocessed_name=self.filenames[self.cur_index]
                row=self.unprocessed_list.row(self.unprocessed_list.findItems(unprocessed_name,Qt.MatchFlag.MatchExactly)[0])
                self.unprocessed_list.takeItem(row)
                
                worker.reset_worker(
                    func,
                    func_params,
                    self.cur_index
                )
                self.cur_index += 1

    def destroy_wokers(self):
        """
        TODO: 下面的代码目前只有一个进程

        需要将其更改为遍历所有的 worker 类，将其全部销毁
        """
        for key,worker in self.workers.items():
            worker.destroy_resource()
            del worker
        self.workers.clear()
        
    def task_over(self):
        self.destroy_wokers()
        self.over_signal.emit()