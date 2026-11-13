import sys

from PySide6.QtWidgets import QApplication

from UI.main_interface import Window

if __name__ == '__main__':
    from multiprocessing import freeze_support
    freeze_support()
    app = QApplication(sys.argv)
    w = Window()
    w.show()
    app.exec()
