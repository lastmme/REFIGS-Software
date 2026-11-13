"""File-centred progress view and non-blocking process scheduler."""
from collections import deque
from datetime import datetime
from multiprocessing import get_context
from pathlib import Path
from queue import Empty
import json
import random
import time
import uuid
import re

import numpy as np
from PySide6.QtCore import QTimer, Signal, Qt, QUrl
from PySide6.QtGui import QDesktopServices, QColor, QTextCursor
from PySide6.QtWidgets import (QWidget, QVBoxLayout, QHBoxLayout, QLabel,
    QPushButton, QProgressBar, QSplitter, QTableWidget, QTableWidgetItem,
    QHeaderView, QPlainTextEdit, QLineEdit, QComboBox, QMessageBox, QSizePolicy)
from pipeline import cal_library, search_args_by_files
from utils.task_runtime import execute_task


class ViewInterface(QWidget):
    over_signal = Signal()

    def __init__(self, object_name, parent=None):
        super().__init__(parent)
        self.setObjectName(object_name)
        self.setStyleSheet('''
            QWidget { font-family: "Microsoft YaHei UI", "Segoe UI"; font-size: 13px; color: #243247; }
            QTableWidget, QPlainTextEdit { background: #ffffff; border: 1px solid #dce4ee; border-radius: 8px; selection-background-color: #e2efff; selection-color: #174c85; }
            QTableWidget { alternate-background-color: #f6f9fd; gridline-color: #edf1f6; }
            QHeaderView::section { background: #f0f5fb; color: #53657d; border: none; padding: 9px 6px; font-weight: 600; }
            QPushButton { background: #ffffff; border: 1px solid #d3ddea; border-radius: 6px; padding: 7px 12px; }
            QPushButton:hover { background: #edf5ff; border-color: #8ab7e8; }
            QPushButton:disabled { color: #98a5b5; background: #f3f5f8; }
            QLineEdit, QComboBox { background: white; border: 1px solid #d3ddea; border-radius: 6px; padding: 7px; }
            QLineEdit:focus { border-color: #3789d4; }
            QProgressBar { background: #e9eff7; border: none; border-radius: 5px; min-height: 14px; text-align: center; font-size: 11px; }
            QProgressBar::chunk { background: #438ed7; border-radius: 5px; }
            QScrollBar:vertical { background: #f1f4f8; width: 10px; margin: 2px; }
            QScrollBar:horizontal { background: #f1f4f8; height: 10px; margin: 2px; }
            QScrollBar::handle { background: #b9c7d7; border-radius: 4px; min-height: 28px; min-width: 28px; }
            QScrollBar::handle:hover { background: #7c9bbd; }
            QScrollBar::add-line, QScrollBar::sub-line { width: 0; height: 0; }
            QScrollBar::add-page, QScrollBar::sub-page { background: transparent; }
            QSplitter::handle { background: #e8eef6; width: 4px; }
        ''')
        self.workers = {}
        self.records = []
        self.state = 'idle'
        self.log_dir = None
        self.selected = None
        self.ctx = get_context('spawn')
        self.timer = QTimer(self)
        self.timer.setInterval(100)
        self.timer.timeout.connect(self._poll)
        root = QVBoxLayout(self)
        root.setContentsMargins(24, 24, 24, 24)
        root.setSpacing(14)
        title = QLabel('分析进度')
        title.setStyleSheet('font-size:24px;font-weight:600')
        root.addWidget(title)
        header = QHBoxLayout()
        self.summary = QLabel('尚未开始任务')
        self.summary.setWordWrap(True)
        self.stop_button = QPushButton('终止任务')
        self.stop_button.setStyleSheet('QPushButton { color: #b42335; } QPushButton:disabled { color: #98a5b5; }')
        self.stop_button.setEnabled(False)
        self.stop_button.clicked.connect(self.stop)
        header.addWidget(self.summary, 1)
        header.addWidget(self.stop_button)
        root.addLayout(header)
        self.progress = QProgressBar()
        self.progress.setValue(0)
        root.addWidget(self.progress)
        self.filter = QComboBox()
        self.filter.addItems(['全部', '待处理', '处理中', '成功', '失败', '已终止', '已取消'])
        self.filter.currentTextChanged.connect(self._filter)
        root.addWidget(self.filter)
        split = QSplitter(Qt.Horizontal)
        self.table = QTableWidget(0, 3)
        self.table.setHorizontalHeaderLabels(['文件 / 阶段', '状态', '耗时'])
        self.table.setAlternatingRowColors(True)
        self.table.setShowGrid(False)
        self.table.setWordWrap(False)
        self.table.setTextElideMode(Qt.ElideMiddle)
        self.table.verticalHeader().hide()
        self.table.verticalHeader().setDefaultSectionSize(40)
        self.table.horizontalHeader().setSectionResizeMode(1, QHeaderView.ResizeToContents)
        self.table.horizontalHeader().setSectionResizeMode(2, QHeaderView.ResizeToContents)
        self.table.horizontalHeader().setSectionResizeMode(0, QHeaderView.Stretch)
        self.table.setSelectionBehavior(QTableWidget.SelectRows)
        self.table.setEditTriggers(QTableWidget.NoEditTriggers)
        self.table.itemSelectionChanged.connect(self._select)
        split.addWidget(self.table)
        pane = QWidget()
        detail = QVBoxLayout(pane)
        self.log_title = QLineEdit('选择文件查看日志')
        self.log_title.setReadOnly(True)
        self.log_title.setMinimumWidth(0)
        self.log_title.setSizePolicy(QSizePolicy.Ignored, QSizePolicy.Fixed)
        self.log_title.setStyleSheet('font-size:15px;font-weight:600;background:#eef4fb;border:none;padding:9px;')
        detail.addWidget(self.log_title)
        tools = QHBoxLayout()
        self.search = QLineEdit()
        self.search.setPlaceholderText('搜索日志，按回车查找下一处')
        self.search.returnPressed.connect(self._find)
        tools.addWidget(self.search)
        for text, callback in [('复制', self._copy), ('到底部', self._bottom), ('打开日志', self._open_log)]:
            button = QPushButton(text)
            button.clicked.connect(callback)
            tools.addWidget(button)
        detail.addLayout(tools)
        self.log = QPlainTextEdit()
        self.log.setReadOnly(True)
        self.log.setMaximumBlockCount(3000)
        self.log.setStyleSheet('font-family:Consolas,"Microsoft YaHei UI";font-size:13px;color:#33465e;background:#fafcff;padding:10px;')
        detail.addWidget(self.log)
        split.addWidget(pane)
        split.setSizes([320, 540])
        root.addWidget(split, 1)
        self.location = QLabel('日志将在结果目录中按运行归档')
        self.location.setWordWrap(True)
        self.location.setStyleSheet('color:#65778e;font-size:12px;')
        self.location.setTextInteractionFlags(Qt.TextSelectableByMouse)
        root.addWidget(self.location)
        folders = QHBoxLayout()
        for text, callback in [('打开日志目录', lambda: self._open(self.records[self.selected]['log'].parent if self.selected is not None else self.log_dir)), ('打开结果目录', lambda: self._open(self.records[self.selected].get('folder', self.result_dir) if self.selected is not None else getattr(self, 'result_dir', None)))]:
            button = QPushButton(text)
            button.clicked.connect(callback)
            folders.addWidget(button)
        folders.addStretch()
        root.addLayout(folders)

    def _open(self, path):
        if path:
            QDesktopServices.openUrl(QUrl.fromLocalFile(str(path)))

    def _open_log(self):
        if self.selected is not None:
            self._open(self.records[self.selected]['log'])

    def _copy(self):
        from PySide6.QtWidgets import QApplication
        QApplication.clipboard().setText(self.log.toPlainText())

    def _bottom(self):
        self.log.moveCursor(QTextCursor.End)

    def _find(self):
        if not self.log.find(self.search.text()):
            self.log.moveCursor(QTextCursor.Start)
            self.log.find(self.search.text())

    def _filter(self):
        for row, record in enumerate(self.records):
            self.table.setRowHidden(row, self.filter.currentText() not in ('全部', record['status']))

    def _select(self):
        row = self.table.currentRow()
        if row < 0 or row >= len(self.records):
            return
        self.selected = row
        record = self.records[row]
        self.log_title.setText(record['name'])
        self.log_title.setCursorPosition(0)
        self.log_title.setToolTip(record['name'])
        try:
            with record['log'].open(encoding='utf-8') as stream:
                self.log.setPlainText(''.join(deque(stream, maxlen=3000)))
        except FileNotFoundError:
            self.log.setPlainText('任务尚未开始')
        except OSError as exc:
            self.log.setPlainText(f'无法读取日志：{exc}')
        self._bottom()

    def _record(self, name, filename):
        row = len(self.records)
        self.records.append(dict(name=name, log=self.log_dir / filename, status='待处理', start=None, end=None))
        self.table.insertRow(row)
        display_name = Path(name).name
        for column, value in enumerate([display_name, '待处理', '—']):
            item = QTableWidgetItem(value)
            item.setToolTip(name)
            self.table.setItem(row, column, item)
        return row

    def _status(self, row, status):
        record = self.records[row]
        record['status'] = status
        if status == '处理中':
            record['start'] = time.monotonic()
        elif record['start'] is not None:
            record['end'] = time.monotonic()
        item = self.table.item(row, 1)
        item.setText(status)
        item.setForeground(QColor({'成功':'#16845b','失败':'#c42b1c','处理中':'#0078d4'}.get(status, '#667085')))
        self._filter()

    def _init_workers(self, function, worker_num, mzxml_files, other_params):
        if self.state in ('running', 'stopping'):
            return
        if not mzxml_files:
            self.over_signal.emit()
            return
        self.run_id = datetime.now().strftime('%Y%m%d_%H%M%S_') + uuid.uuid4().hex[:8]
        self.result_dir = Path(other_params['result_dir']) / self.run_id
        self.log_dir = self.result_dir / '_shared' / 'logs'
        task_params = other_params.copy()
        task_params['result_dir'] = str(self.result_dir / '_shared')
        try:
            self.log_dir.mkdir(parents=True)
            (self.log_dir / 'run.log').write_text(json.dumps(dict(run_id=self.run_id, files=mzxml_files, params=task_params), ensure_ascii=False, indent=2), encoding='utf-8')
        except OSError as exc:
            QMessageBox.critical(self, '无法创建日志', str(exc))
            self.over_signal.emit()
            return
        self.records = []
        self.table.setRowCount(0)
        self.selected = None
        self.log.clear()
        self.function, self.params = function, task_params
        self.files = list(mzxml_files)
        self.worker_num = max(1, min(worker_num, len(self.files)))
        self.pending = deque(range(len(self.files)))
        self.cancel = self.ctx.Event()
        self.state = 'running'
        self.started = time.monotonic()
        self._record('谱图库准备', 'library.log')
        self.auto = self.params['is_auto_search_refigs_params']
        if self.auto:
            self._record('自动参数搜索', 'parameter_search.log')
        self.offset = len(self.records)
        for index, path in enumerate(self.files):
            self._record(str(Path(path)), f'{index:04d}.log')
            safe_name = re.sub(r'[^\w.\-]+', '_', Path(path).stem)[:80]
            folder = self.result_dir / f'{index:04d}_{safe_name}'
            try:
                folder.mkdir()
            except OSError as exc:
                self._fail_batch(str(exc))
                self._finish()
                return
            self.records[-1]['folder'] = folder
            self.records[-1]['log'] = folder / 'task.log'
        self.location.setText(f'任务目录：{self.result_dir}\n日志目录：{self.log_dir}')
        self.stop_button.setEnabled(True)
        self.progress.setRange(0, 0)
        self._launch(0, cal_library, self.params.copy(), 'library')
        self.table.selectRow(0)
        self.timer.start()

    def _launch(self, row, function, params, phase):
        queue = self.ctx.Queue()
        process = self.ctx.Process(target=execute_task, args=(function, params, str(self.records[row]['log']), queue, self.cancel))
        try:
            process.start()
        except Exception as exc:
            queue.close()
            self._status(row, '失败')
            self._fail_batch(str(exc))
            return
        self.workers[row] = dict(process=process, queue=queue, phase=phase, result=None, forced=False)
        self._status(row, '处理中')

    def _dispatch(self):
        while self.state == 'running' and self.pending and len(self.workers) < self.worker_num:
            index = self.pending.popleft()
            params = self.params.copy()
            params.pop('precomputed_library_dir', None)
            params['mzxml_files'] = [self.files[index]]
            params['library_dir'] = self.params['result_dir']
            params['result_dir'] = str(self.records[index + self.offset]['folder'])
            self._launch(index + self.offset, self.function, params, 'analysis')

    def _poll(self):
        for row, worker in list(self.workers.items()):
            drained = False
            for _ in range(200):
                try:
                    event = worker['queue'].get_nowait()
                except Empty:
                    drained = True
                    break
                if event['type'] == 'result':
                    worker['result'] = event
                elif self.selected == row:
                    bar = self.log.verticalScrollBar()
                    follow = bar.value() >= bar.maximum() - 4
                    cursor = self.log.textCursor()
                    cursor.movePosition(QTextCursor.End)
                    cursor.insertText(event['text'])
                    if follow:
                        self._bottom()
            process = worker['process']
            if process.is_alive():
                if self.state == 'stopping' and time.monotonic() - self.stop_time > 3:
                    process.kill()
                    worker['forced'] = True
                continue
            if not drained:
                continue
            # A joined process has flushed its queue feeder; drain remaining events
            # before interpreting an absent terminal event on the next tick.
            if not worker.get('exited'):
                process.join(timeout=0)
                worker['exited'] = True
                continue
            result = worker['result'] or dict(status='failed', detail=f'子进程异常退出，退出码 {process.exitcode}')
            worker['queue'].close()
            worker['queue'].cancel_join_thread()
            process.close()
            del self.workers[row]
            if self.state == 'stopping':
                self._status(row, '已终止')
            else:
                self._status(row, '成功' if result['status'] == 'success' else '失败')
            try:
                with self.records[row]['log'].open('a', encoding='utf-8') as stream:
                    stream.write(f'\n[调度器] {self.records[row]["status"]}，强制终止={worker["forced"]}\n{result["detail"]}\n')
            except OSError as exc:
                self._fail_batch(f'日志写入失败：{exc}')
            if self.selected == row:
                self._select()
            if self.state != 'running':
                continue
            if worker['phase'] != 'analysis' and result['status'] != 'success':
                self._fail_batch(result['detail'])
            elif worker['phase'] == 'library' and self.auto:
                params = self.params.copy()
                params['mzxml_files'] = random.Random(params.get('seed', 13)).sample(self.files, min(5, len(self.files)))
                self._launch(1, search_args_by_files, params, 'search')
            else:
                if worker['phase'] == 'search':
                    try:
                        self.params.update(np.load(Path(self.params['result_dir']) / 'refigs_params.npy', allow_pickle=True).item())
                    except Exception as exc:
                        self._fail_batch(str(exc))
                        continue
                self.progress.setRange(0, len(self.files))
                self._dispatch()
        self._refresh()
        if not self.workers and self.state in ('running', 'stopping'):
            self._finish()

    def _refresh(self):
        counts = {status: sum(r['status'] == status for r in self.records[self.offset:]) for status in ['成功','失败','已终止','已取消']}
        ended = sum(counts.values())
        if self.progress.maximum():
            self.progress.setValue(ended)
        phase = next(({'library':'准备谱图库', 'search':'自动寻参', 'analysis':'分析文件'}[w['phase']] for w in self.workers.values()), '')
        self.summary.setText(f'{"正在终止" if self.state == "stopping" else "运行中"} {phase} · 已结束 {ended}/{len(self.files)} · 成功 {counts["成功"]} · 失败 {counts["失败"]} · {int(time.monotonic()-self.started)} 秒')
        for row, record in enumerate(self.records):
            if record['start'] is not None:
                seconds = int((record['end'] or time.monotonic()) - record['start'])
                self.table.item(row, 2).setText(f'{seconds} 秒')

    def stop(self):
        if self.state != 'running':
            return
        self.state = 'stopping'
        self.stop_time = time.monotonic()
        self.cancel.set()
        self.pending.clear()
        self.stop_button.setEnabled(False)
        for row, record in enumerate(self.records):
            if record['status'] == '待处理':
                self._status(row, '已取消')

    def _fail_batch(self, detail):
        self.stop()
        self.batch_error = detail
        self.location.setText(f'任务失败：{detail}\n日志：{self.log_dir}')

    def _finish(self):
        self.timer.stop()
        failed = any(r['status'] == '失败' for r in self.records)
        self.state = 'failed' if failed else ('cancelled' if self.state == 'stopping' else 'complete')
        self.summary.setText(self.summary.text().replace('正在终止', '已终止').replace('运行中', '失败' if failed else '已完成'))
        self.stop_button.setEnabled(False)
        self.progress.setRange(0, len(self.files))
        self.progress.setValue(sum(r['status'] != '待处理' for r in self.records[self.offset:]))
        self.over_signal.emit()

    def destroy_wokers(self):
        # Compatibility entry point: window close now waits asynchronously.
        self.stop()
