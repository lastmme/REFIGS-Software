from pathlib import Path
import json
import multiprocessing
import time


def calculation(output):
    import numpy as np
    from scipy.optimize import nnls
    import cvxopt
    from numba import njit
    from sklearn.discriminant_analysis import LinearDiscriminantAnalysis
    import pandas
    from lxml import etree
    from pyteomics import mzxml
    import statsmodels.api
    solve = njit(lambda x: x.sum())
    assert solve(np.array([1., 2.])) == 3.
    assert nnls(np.eye(2), np.array([1., 2.]))[1] == 0
    Path(output).write_text('worker OK', encoding='utf-8')


def run(output):
    import os
    os.environ['QT_QPA_PLATFORM'] = 'offscreen'
    from PySide6.QtWidgets import QApplication
    from UI.main_interface import Window
    app = QApplication([])
    window = Window()
    worker_output = str(Path(output).with_suffix('.worker.txt'))
    process = multiprocessing.get_context('spawn').Process(target=calculation, args=(worker_output,))
    process.start()
    deadline = time.monotonic() + 120
    while process.is_alive() and time.monotonic() < deadline:
        app.processEvents()
        time.sleep(.05)
    if process.is_alive():
        process.kill()
        process.join(5)
        raise RuntimeError('frozen worker timeout')
    assert process.exitcode == 0, f'worker exit={process.exitcode}'
    assert Path(worker_output).read_text(encoding='utf-8') == 'worker OK'
    window.close()
    Path(output).write_text(json.dumps({'ui': 'OK', 'spawn': 'OK', 'scientific_libraries': 'OK'}), encoding='utf-8')
