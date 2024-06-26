# -*- coding: utf-8 -*-
"""
Created on Tue Jan 16 16:44:39 2024

@author: Lab Di Martino
"""


"""
Created on Wed Sep 13 16:35:28 2023

@author: Lab Di Martino
Trung - xtn20
January 2024
Interface and control for Keithley


"""


from PyQt5 import QtWidgets, uic
from PyQt5.QtCore import QTimer, QTime, Qt, QRunnable, QThreadPool
from nplab.instrument.camera.Andor.andor_sdk_rig2 import AndorBase
from nplab.instrument.spectrometer.Kymera import Kymera
from nplab.experiment import Experiment, ExperimentStopped

import sys
import os
import numpy as np
import matplotlib.pyplot as plt
import threading
import matplotlib.colors as colors
from nplab.utils.gui import QtCore, QtGui, uic, get_qt_app, show_widget
from nplab.ui.ui_tools import UiTools
import ctypes
from ctypes import byref, c_int, c_ulong, c_double
import pyqtgraph as pg
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg
from matplotlib.figure import Figure
import matplotlib.pyplot as plt
from datetime import datetime
import time

class testandor(QtWidgets.QWidget, UiTools):
    #inherit QWidget methods
    def __init__(self, ui_file = os.path.join(os.path.dirname(__file__),'xyz_stack_DF_Raman.ui')):
        super(testandor, self).__init__() 
        uic.loadUi(ui_file, self)
        
app = QtWidgets.QApplication(sys.argv)
main = testandor()
main.show()
sys.exit(app.exec_())

        












 



































        
        