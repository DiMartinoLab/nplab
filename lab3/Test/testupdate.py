# -*- coding: utf-8 -*-
"""
Created on Wed Nov  1 15:59:34 2023

@author: Lab Di Martino
"""

import nplab
from PyQt5 import QtWidgets, uic
from PyQt5.QtCore import QTimer, QTime, Qt, QRunnable, QThreadPool
import sys
#import os
import numpy as np

from nplab.ui.ui_tools import UiTools
import time




class OlympusCamera(QtWidgets.QWidget, UiTools):
    #inherit QWidget methods
    def __init__(self, ui_file = 'testupdate.ui',   parent=None):
        
        super(OlympusCamera, self).__init__() 
        uic.loadUi(ui_file, self)
#        self.doubleSpinBox.valueChanged.connect(self.update_text)
        self.pushButton.clicked.connect(self.Capture)
        self.threadCount = QThreadPool.globalInstance().maxThreadCount()
        self.pool = QThreadPool.globalInstance()
        
        self.YLim1_box.valueChanged.connect(self.update_figure_lim)
        self.YLim2_box.valueChanged.connect(self.update_figure_lim)
        self.XLim1_box.valueChanged.connect(self.update_figure_lim)
        self.XLim2_box.valueChanged.connect(self.update_figure_lim)
        
        self.wa =   [70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80]
        self.spec = [1, 5, 8, 2, 1, 7, 8, 4, 2, 6, 1]
        self.Spectra_display.canvas.ax.plot(self.wa, self.spec)
        self.update_figure_lim()
        self.Spectra_display.canvas.draw()
        
    def update_text(self):
        self.lineEdit.setText(str(self.doubleSpinBox.value()))
        print('replied')

    def Capture(self):
        start_capture = sub_capture(self.Spectra_display)
        self.pool.start(start_capture)
        
    def update_figure_lim(self):
        self.Spectra_display.minY = min(float(self.YLim1_box.value()), float(self.YLim2_box.value()))
        self.Spectra_display.maxY = max(float(self.YLim1_box.value()), float(self.YLim2_box.value()))
        if self.Spectra_display.minY == self.Spectra_display.maxY:
            self.Spectra_display.minY = min(self.spec)
            self.Spectra_display.maxY = max(self.spec)   
        self.Spectra_display.canvas.ax.set_ylim([self.Spectra_display.minY, self.Spectra_display.maxY])
        
        """ Display WN axis"""
        self.Spectra_display.minX = min(float(self.XLim1_box.value()), float(self.XLim2_box.value()))
        self.Spectra_display.maxX = max(float(self.XLim1_box.value()), float(self.XLim2_box.value()))
        self.Spectra_display.canvas.ax.set_xlim([self.Spectra_display.minX, self.Spectra_display.maxX])
        

        self.Spectra_display.canvas.ax2.set_xlim([self.Spectra_display.minX, self.Spectra_display.maxX])
#        
        self.laser_wl = 74
        
        def WL2WN(wavelength):
            return 1e4/(self.laser_wl) - 1e4/(wavelength)
        
        def WN2WL(k):
            return 1e4/(1e4/(self.laser_wl) - k)
        
        kmin = WL2WN(self.Spectra_display.minX)
        kmax = WL2WN(self.Spectra_display.maxX)
#        print(kmin)
#        print(kmax)
        scale = 1
        kticksn = list(range(0, round(kmin)*scale, -2))
        kticksn.reverse()
        kticksp = list(range(0, round(kmax)*scale, 2))
        kticks = [*kticksn, *kticksp]
        print(kticks)
        kticksInNm = WN2WL(np.array(kticks)/scale)
        print(kticksInNm)
        kticksstr = []
#        print(self.wa)
#        print(kticks)
        for i in kticks:
            kticksstr.append(str(i/scale))
        self.Spectra_display.canvas.ax2.set_xticks(kticksInNm)
        self.Spectra_display.canvas.ax2.set_xticklabels(kticksstr)
        self.Spectra_display.canvas.ax2.set_xlim([self.Spectra_display.minX, self.Spectra_display.maxX])
#        """ Display WN axis"""
#        
#        for label in (self.Spectra_display.canvas.ax.get_xticklabels() +  \
#                      self.Spectra_display.canvas.ax.get_yticklabels() + \
#                      self.Spectra_display.canvas.ax2.get_xticklabels()):
#            label.set_fontsize(16)
        self.Spectra_display.canvas.draw()
        
class sub_capture(QRunnable):
    def __init__(self, Spectra_display):
        super().__init__()
        i = 0
        self.Spectra_display = Spectra_display
        
    def run(self):
        wa =   [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
        spec = [1, 5, 8, 2, 1, 7, 8, 4, 2, 6, 1]
        spec = [i * 5 for i in spec]
        self.Spectra_display.canvas.ax.cla()
        self.Spectra_display.canvas.ax.plot(wa, spec)

        if self.Spectra_display.minY == self.Spectra_display.maxY:
            self.Spectra_display.minY = min(spec)
            self.Spectra_display.maxY = max(spec)   
            
        self.Spectra_display.canvas.ax.set_ylim([self.Spectra_display.minY, self.Spectra_display.maxY])
        self.Spectra_display.canvas.ax.set_xlim([self.Spectra_display.minX, self.Spectra_display.maxX])
        self.Spectra_display.canvas.ax2.set_xlim([self.Spectra_display.minX, self.Spectra_display.maxX])
        self.Spectra_display.canvas.draw()
        
app = QtWidgets.QApplication(sys.argv)
main = OlympusCamera()
main.show()
sys.exit(app.exec_())