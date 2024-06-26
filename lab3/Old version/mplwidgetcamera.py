# -*- coding: utf-8 -*-
"""
Created on Fri Nov 18 13:19:46 2022

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
        
#        self.rcParams['axes.linewidth'] = 1
#        self.rcParams['xtick.major.width'] = 0
#        self.rcParams['ytick.major.width'] = 0
#        self.rcParams['xtick.major.size'] = 0
#        self.rcParams['ytick.major.size'] = 0
#        self.rcParams['xtick.top'] = 'False'
#        self.rcParams['ytick.right'] = 'False'
#        self.rcParams['font.size'] = '12'
#        self.tick_params(which = 'both', direction = "in")
        # self.fig.auto_layout()
        Canvas.__init__(self, self.fig)
        Canvas.setSizePolicy(self, QtWidgets.QSizePolicy.Preferred, QtWidgets.QSizePolicy.Preferred)
        Canvas.updateGeometry(self)

# Matplotlib widget
class MplWidgetCamera(QtWidgets.QWidget):
    def __init__(self, parent=None):
        QtWidgets.QWidget.__init__(self, parent)   # Inherit from QWidget
        self.canvas = MplCanvas()                  # Create canvas object
        self.vbl = QtWidgets.QVBoxLayout()         # Set box for plotting
        self.vbl.addWidget(self.canvas)
        self.setLayout(self.vbl)
        

    
    
    