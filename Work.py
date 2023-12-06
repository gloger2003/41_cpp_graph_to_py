import math
from typing import List
from PyQt6.QtWidgets import *
from PyQt6.QtGui import *
from PyQt6.QtCore import *

import numpy as np

class MyDecision:
    alpha: list[float]
    betta: list[float]

    psi: List[float]

    mua: float
    va: float
    mub: float
    vb: float
    
    k1: float
    k2: float
    k3: float
    k4: float
    
    y0: float
    stepX: float
    x0: float
    xl: float
    yl: float
    new_y0: float
    new_yl: float
    new_x0: float
    new_xl: float
        
    delta: float 
    Ld: float
    F0: float
    nd: float
    N0: float
    
    points: int
    
    def __init__(self, step: float, my_y0: float, my_x0: float,
                 my_yl: float, my_xl: float, e: float, Nd: float,
                 _delta: float, num: int):
        """
        my_decision(double step, double my_y0, double my_x0, double my_yl, double my_xl, double e, double Nd, double _delta, int num)
        """
        self.y0 = my_y0
        self.x0 = my_x0
        self.xl = my_xl
        self.yl = my_yl
        self.delta = _delta
        self.points = num
        
        self.T: float = 300
        self.kb: float = 1.38064852e-23
        self.qe: float = 1.6e-19
        self.e0: float = 8.86e-12
        
        self.F0 = (self.kb * self.T) / self.qe
        self.new_y0 = self.y0 / self.F0
        self.new_yl = self.yl / self.F0

        self._Nd: float = step
        self.nd = Nd
        self.N0 = self._Nd / self.nd
        self.Ld = math.sqrt((e * self.e0 * self.kb * self.T) / (self.qe * self.qe * self._Nd))
        self.new_x0 = self.x0 / self.Ld
        self.new_xl = self.xl / self.Ld
        self.stepX = (self.new_xl) / self.points

        # vector<double>(points + 1)
        # FIXME
        self.psi: List[float] = []

    def funcAlpha(self, x: float, y: float, _y0: float) -> float:
        """ double funcAlpha(double x, double y, double _y0); """
        return y * y * self.funcQ(_y0) + 1

    def funcBetta(self, x: float, y: float, alpha: float, _y0: float) -> float:
        """ double funcBetta(double x, double y, double alpha, double _y0); """
        return -alpha * self.funcR(_y0) + alpha * y * self.funcQ(_y0)

    def progonka(self, x: List[float], y: List[float], vecy0: List[float]) -> None:
        """ void progonka(std::vector<double>& x, std::vector<double>& y, std::vector<double> vecy0); """
        self.Eyler(x, y, 0, vecy0)
        self.Eyler(x, y, 1, vecy0)

    def control(self, x: List[float], y: List[float]) -> None:
        """ void control(std::vector<double>& x, std::vector<double>& y); """
        x.clear()
        
        vecy0: List[float]
        my_try: List[float]
        norma: float
        iter_: int = 0
        
        self.SetY0(y, x)
        my_try = y

        def do_iter():
            global norma
            global iter_
            
            vecy0 = y
            
            self.progonka(x, y, vecy0)
            
            norma = self.Norma(y, vecy0)
            
            iter_ += 1
            if iter_ > 100:
                # MessageBox(NULL, L"iter > 100", L"iter", NULL);
                return False
            return True
        
        if do_iter():
            while norma > self.delta:
                res = do_iter()
                if not res:
                    break
                
        # y = alpha;
        self.Y_obr(x, y)

    def SetY0(self, y: List[float], x: List[float]) -> None:
        """ void SetY0(std::vector<double>& y, std::vector<double>& x); """
        self.new_y0 = self.y0
        self.new_yl = self.yl
        y.append(self.new_y0)
        x.append(self.new_x0)
        
        for i in range(1, self.points + 1, 1):
            y.append(self.new_y0 + (self.new_yl - self.new_y0) * i / (self.points));
            x.append(self.new_xl * i / self.points);
        pass

    def Norma(self, y: List[float], y0: List[float]) -> float:
        """ double Norma(std::vector<double> y, std::vector<double> y0); """
        razn: float = 0

        if y.size() != self.vecy0.size():
            # FIXME
            # MessageBox(NULL, L"Error", L"error y y0", NULL);
            return 0.0

        for i in range(self.points + 1):
            razn += (y[i] - self.vecy0[i]) * (y[i] - self.vecy0[i])
            
        return razn

    def Y_obr(self, x: List[float], y: List[float]) -> None:
        """ void Y_obr(std::vector<double>& x, std::vector<double>& y); """
        if len(y) != len(x):
            # FIXME
            # MessageBox(NULL, L"Error", L"error x y", NULL);
            return

        for i in range(self.points + 1):
            x[i] *= self.Ld * 1e6
            # //y[i] *= F0;
        pass

    def func(self, x: float, y: float, vy: float, my_y0: float) -> float:
        """ double func(double x, double y, double vy, double my_y0); """
        return self.nd * (1. - math.exp(-y))

    def funcQ(self, y0: float) -> float:
        """ double funcQ(double y0); """
        return -self.nd * math.exp(-y0)

    def funcR(self, y0: float) -> float:
        """ double funcR(double y0); """
        return self.func(y0, y0, y0, y0) + self.funcQ(y0) * y0

    def Eyler(self, x: List[float], vec1: List[float], function: int, vecy0: List[float]) -> None:
        """ void Eyler(std::vector<double>& x, std::vector<double>& vec1, int function, std::vector<double> vecy0); """
        y: List[float]

        if function == 0:
            self.alpha = [0] * (self.points + 1)
            self.betta = [0] * (self.points + 1)

            self.alpha[0] = 0
            self.betta[0] = self.y0

            # for (int i = 1; i <= points; i++):
            for i in range(1, self.points, 1):
                self.alpha[i] = self.alpha[i - 1] + (self.alpha[i - 1]
                                                     * self.alpha[i - 1]
                                                     * self.funcQ(vecy0[i - 1]) + 1) * self.stepX
                self.betta[i] = self.betta[i - 1] + (self.alpha[i - 1]
                                                     * self.betta[i - 1]
                                                     * self.funcQ(vecy0[i - 1])
                                                     - self.alpha[i - 1]
                                                     * self.funcR(vecy0[i - 1])) * self.stepX
                
            vec1[self.points] = self.new_yl
            self.psi[self.points] = (vec1[self.points] - self.betta[self.points]) / self.alpha[self.points]
        else:
            for i in range(self.points - 1, 0, -1):
                self.psi[i] = self.psi[i + 1] - (self.funcR(vecy0[i + 1])
                                                 - self.funcQ(vecy0[i + 1])
                                                 * (self.alpha[i + 1] 
                                                    * self.psi[i + 1]
                                                    + self.betta[i + 1])) * self.stepX

            # FIXME
            # vec1 = vector<double>(points + 1);
            vec1.clear()
            for i in range(self.points):
                vec1.append(0)
                vec1[i] = self.alpha[i] * self.psi[i] + self.betta[i]