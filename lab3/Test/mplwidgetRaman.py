# -*- coding: utf-8 -*-
"""
Created on Wed Sep 13 12:58:13 2023

@author: Trung
"""

# Imports
from PyQt5 import QtWidgets
from matplotlib.pyplot import Figure
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as Canvas
import matplotlib

# Ensure using PyQt5 backend
matplotlib.use('QT5Agg')

# Matplotlib canvas class to create figure
class MplCanvas(Canvas):
    def __init__(self):
        self.fig = Figure()
        self.ax = self.fig.add_subplot(111)
        self.ax.set_ylabel('y-label')
        self.ax.set_xlabel('x-label')
        self.ax2 = self.ax.twiny()
        self.ax.set_ylabel('y-label')
        self.ax.tick_params(axis = "y", direction = "in")
        self.ax.tick_params(axis = "x", direction = "in")
        self.ax2.tick_params(axis = "x", direction = "in", top = 'true')
        
        # self.fig.auto_layout()
        Canvas.__init__(self, self.fig)
        Canvas.setSizePolicy(self, QtWidgets.QSizePolicy.Preferred, QtWidgets.QSizePolicy.Preferred)
        Canvas.updateGeometry(self)

# Matplotlib widget mplwidgetRaman
class mplwidgetRaman(QtWidgets.QWidget):
    def __init__(self, parent=None):
        QtWidgets.QWidget.__init__(self, parent)   # Inherit from QWidget
        self.canvas = MplCanvas()                  # Create canvas object
        self.vbl = QtWidgets.QVBoxLayout()         # Set box for plotting
        self.vbl.addWidget(self.canvas)
        self.setLayout(self.vbl)
        
class pltsetting:
    def setting(self):
        self.rcParams['axes.linewidth'] = 1
        self.rcParams['xtick.major.width'] = 1
        self.rcParams['ytick.major.width'] = 1
        self.rcParams['xtick.major.size'] = 5
        self.rcParams['ytick.major.size'] = 5
        self.rcParams['xtick.top'] = 'True'
        self.rcParams['ytick.right'] = 'True'
        self.rcParams['font.size'] = '16'
        self.tick_params(which = 'both', direction = "in")
        # self.set_linewidth = 5
        
    def plotsmt(self, x_range, x_label, y_range, y_label, legendloc):
        """
        x/y_range = [bottom, upper], 0 for no limitation \n
        x/y_label = x/y-axis label \n
        legendloc = legend location, '' for no legend
        """
        if len(x_range)>1:
            self.xlim(x_range)
        self.xlabel(x_label)
        if len(y_range)>1:
            self.ylim(y_range)
        self.ylabel(y_label)
        if len(legendloc)>0:
            self.legend(loc=legendloc, fontsize=10)
        self.show()
    
    def plotticks(self, xticks, yticks):
        """ x/yticks = [min, max, increasment], 0 for no limitation """
        import numpy as np
        if len(xticks)>1:
            self.xticks(np.arange(xticks[0], xticks[1], xticks[2]))
        if len(yticks)>1:
            self.yticks(np.arange(yticks[0], yticks[1], yticks[2]))
    
    
    