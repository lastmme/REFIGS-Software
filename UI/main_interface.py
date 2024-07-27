import sys

import pixmaps
from PySide6.QtWidgets import QFrame, QHBoxLayout, QApplication
from PySide6.QtCore import Qt, Signal
from PySide6.QtGui import QIcon
from qfluentwidgets import (
    NavigationItemPosition,
    MSFluentWindow,
    SubtitleLabel,
    setFont
)

from pipeline import pipeline
from .view_interface import ViewInterface
from .software_introduction import IntroducationInterface
from .setting_interface import SettingInterface


class Widget(QFrame):
    def __init__(self, text: str, parent=None):
        super().__init__(parent=parent)
        self.label = SubtitleLabel(text, self)
        self.hBoxLayout = QHBoxLayout(self)

        setFont(self.label, 24)
        self.label.setAlignment(Qt.AlignmentFlag.AlignCenter)
        self.hBoxLayout.addWidget(self.label, 1, Qt.AlignmentFlag.AlignCenter)
        self.setObjectName(text.replace(' ', '-'))


class Window(MSFluentWindow):

    def __init__(self):
        super().__init__()

        self.setting_interface = SettingInterface('Setting Interface', self)
        self.view_progress_interface = ViewInterface('View Progress', self)
        self.introduction_interface = IntroducationInterface(
            'Software Introduction', self)
        self.init_window()
        self.init_navigation()
        self._connect_signal2solt()

    def _connect_signal2solt(self):
        self.setting_interface.excute_signal.connect(self.run)
        self.view_progress_interface.over_signal.connect(self.task_over)

    def run(self):
        """
        TODO:
        这个地方 setworker 可能需要该成 init_worker 函数

        相当于初始化 worker，因为后面还需要进行 worker 的调度
        """
        params = self.setting_interface.get_params()
        worker_num=params[0]
        mzxml_files=params[1]['mzxml_files']
        other_params=params[1]
        del other_params['mzxml_files']
        self.view_progress_interface._init_workers(pipeline,worker_num,mzxml_files,other_params)
        self.navigationInterface.setCurrentItem(
            self.setting_interface.objectName()
        )
        self.switchTo(self.view_progress_interface)
        # self.view_progress_interface.callib_first()
        
    def task_over(self):
        self.setting_interface.reset()

    def init_navigation(self):
        self.addSubInterface(
            self.setting_interface,
            ':/normal/svgs/normal/setting.svg',
            '参数配置',
            ':/selected/svgs/selected/setting.svg'
        )
        self.addSubInterface(
            self.view_progress_interface,
            ':/normal/svgs/normal/view_progress.svg',
            '查看进度',
            ':/selected/svgs/selected/view_progress.svg'
        )
        self.addSubInterface(
            self.introduction_interface,
            ':/normal/svgs/normal/Introduction.svg',
            '软件简介',
            ':/selected/svgs/selected/Introduction.svg',
            NavigationItemPosition.BOTTOM
        )

    def init_window(self):
        self.resize(900, 700)
        self.setWindowIcon(QIcon(':/selected/svgs/RE-FIGS.svg'))
        self.setWindowTitle('RE-FIGS')

        desktop = QApplication.screens()[0].availableGeometry()
        w, h = desktop.width(), desktop.height()
        self.move(w//2 - self.width()//2, h//2 - self.height()//2)

    def closeEvent(self, e):
        self.view_progress_interface.destroy_wokers()
        return super().closeEvent(e)


if __name__ == '__main__':
    app = QApplication(sys.argv)
    w = Window()
    w.show()
    app.exec()
