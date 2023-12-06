import sys
from PyQt6.QtWidgets import *
from PyQt6.QtGui import *
from PyQt6.QtCore import *

import MyGraph

import pyqtgraph as pg


class MainWindow(QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle('Legirovanie')

        self.setCentralWidget(MyGraph.MyGraph(self))
        
        
if __name__ == '__main__':
    # 1275x530
    app = QApplication([])
    
    window = MainWindow()
    window.setFixedSize(1275, 530)
    window.show()
    
    sys.exit(app.exec())