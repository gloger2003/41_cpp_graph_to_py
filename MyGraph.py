from pprint import pprint
from PyQt6.QtWidgets import *
from PyQt6.QtGui import *
from PyQt6.QtCore import *

import pyqtgraph as pg
import numpy as np
from Work import MyDecision

import main


class MyGraph(QWidget):
    def __init__(self, window: 'main.MainWindow'):
        super().__init__(window)
        
        self.plot_graph = pg.PlotWidget()
        self.plot_graph.setBackground("w")

        self.main_vbox = QVBoxLayout(self)
        self.main_vbox.addWidget(self.plot_graph)
        
        step = 0.1
        
        # TODO: From UI
        Nd = 1e21
        e = 12
        V = 2
        delta = 1e-8
        number_point = 500
        nd = 1e-1
        L = 2
        
        examp = MyDecision(Nd, V, 0, 0, L * 1e-6, e, nd, delta, number_point)
        
        # x: np.ndarray = np.ndarray((number_point), dtype=np.double)
        # y: np.ndarray = np.zeros((number_point), dtype=np.double)
        
        x = []
        y = []
        
        examp.control(x, y)
        
        pen = pg.mkPen(color=(255, 0, 0))
        # time = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
        # temperature = [30, 32, 34, 32, 33, 31, 29, 32, 35, 30]
        yL = [L / 13 * i for i in range(L)]
        
        xL = [L / 100 * i / (number_point / 100) for i in range(len(x))][1:]
        # pprint(yn)
        # yn = []
        # xn = []
        # for k in range(len(xL)):
        #     yn.append(y[k])
        #     xn.append(x[k])
        #     if xL[k] == L:
        #         return
            
        x = [L / 100 * i / (number_point / 100) for i in range(len(x))]
        self.plot_graph.plot(x, y, pen=pen)