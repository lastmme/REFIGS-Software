import sys
import os
from pathlib import Path

if getattr(sys, 'frozen', False):
    log_dir = Path(os.environ.get('LOCALAPPDATA', str(Path.home()))) / 'RE-FIGS'
    log_dir.mkdir(parents=True, exist_ok=True)
    if sys.stdout is None:
        sys.stdout = open(os.devnull, 'w')
    if sys.stderr is None:
        sys.stderr = open(log_dir / 'startup.log', 'a', encoding='utf-8', buffering=1)

from PySide6.QtWidgets import QApplication

from UI.main_interface import Window

if __name__ == '__main__':
    from multiprocessing import freeze_support
    freeze_support()
    if '--smoke-test' in sys.argv:
        from packaging_smoke import run
        run(sys.argv[sys.argv.index('--smoke-test') + 1])
        sys.exit(0)
    app = QApplication(sys.argv)
    w = Window()
    w.show()
    app.exec()
