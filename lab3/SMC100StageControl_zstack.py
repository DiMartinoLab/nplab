# -*- coding: utf-8 -*-
"""
Created on Wed Feb 15 10:39:38 2023

@author: Lab Di Martino
"""

import nplab
#from nplab.instrument.spectrometer.kandor import Kandor
#from nplab.instrument.spectrometer.seabreeze import OceanOpticsSpectrometer, OceanOpticsControlUI
#from nplab.instrument.electronics.keithley_2636b_smu import Keithley2636B as Keithley
# from nplab.instrument.stage.smaract_mcs import SmaractMCSSerial
#from nplab.instrument.shutter.Arduino_ttl_shutter import Arduino_tri_shutter as shutter
from nplab.instrument.light_sources.matchbox_laser import MatchboxLaser
from nplab.instrument.stage.SMC100_lib_zstack import SMC100
from PyQt5.QtCore import QRunnable, Qt, QThreadPool
import nplab.utils.gui 
import nplab.datafile as datafile
#from nplab.instrument.spectrometer import Spectrometer
from nplab.experiment import Experiment, ExperimentStopped
from nplab.utils.notified_property import DumbNotifiedProperty, NotifiedProperty
from nplab.utils.gui import QtCore, QtGui, uic, get_qt_app, show_widget, popup_widget
from nplab.ui.ui_tools import UiTools
import os
import time
import threading
import numpy as np
import pyqtgraph as pg
from pyqtgraph import PlotWidget, plot
from PyQt5 import QtWidgets, uic

from nplab.ui.ui_tools import QtWidgets
import sys
# never wait for more than this e.g. during wait_states
MAX_WAIT_TIME_SEC = 300

# time to wait after sending a command. This number has been arrived at by
# trial and error
COMMAND_WAIT_TIME_SEC = 1

# States from page 65 of the manual
STATE_NOT_REFERENCED_FROM_RESET = '0A'
STATE_NOT_REFERENCED_FROM_CONFIGURATION = '0C'
STATE_READY_FROM_HOMING = '32'
STATE_READY_FROM_MOVING = '33'

STATE_CONFIGURATION = '14'

STATE_DISABLE_FROM_READY = '3C'
STATE_DISABLE_FROM_MOVING = '3D'
STATE_DISABLE_FROM_JOGGING = '3E'


class SMC100ReadTimeOutException(Exception):
    def __init__(self):
        super(SMC100ReadTimeOutException, self).__init__('Read timed out')


class SMC100WaitTimedOutException(Exception):
    def __init__(self):
        super(SMC100WaitTimedOutException, self).__init__('Wait timed out')


class SMC100DisabledStateException(Exception):
    def __init__(self, state):
        super(SMC100DisabledStateException, self).__init__('Disabled state encountered: ' + state)


class SMC100RS232CorruptionException(Exception):
    def __init__(self, c):
        super(SMC100RS232CorruptionException, self).__init__('RS232 corruption detected: %s' % (hex(ord(c))))


class SMC100InvalidResponseException(Exception):
    def __init__(self, cmd, resp):
        s = 'Invalid response to %s: %s' % (cmd, resp)
        super(SMC100InvalidResponseException, self).__init__(s)

""" BODY """
class SMC100_window_object(QtWidgets.QWidget, UiTools):
    #inherit QWidget methods
    """ use this line for lab3 experiment"""
    def __init__(self, SMC100_instance, ui_file = os.path.join(os.path.dirname(__file__),'SMC100StageControl_zstack.ui'),  parent=None):
#    """ use this line for running SMC100 only"""
#    def __init__(self, ui_file = os.path.join(os.path.dirname(__file__),'SMC100StageControl_zstack.ui'),  parent=None):
        super(SMC100_window_object, self).__init__() 
        uic.loadUi(ui_file, self)
        """ use this line for lab3 experiment"""    
        self.SMC100 = SMC100_instance
        """ use this line for running SMC100 only"""
#        self.SMC100 = SMC100('COM1', (1,2,3))
        
        self.InitialiseButton.clicked.connect(self.initialise_stage)
        self.DisconnectButton.clicked.connect(self.disconnect_stage)
        
        self.stage_mup_button.clicked.connect(self.SM100_move_up)
        self.stage_mdown_button.clicked.connect(self.SM100_move_down)
        self.stage_mleft_button.clicked.connect(self.SM100_move_left)
        self.stage_mright_button.clicked.connect(self.SM100_move_right)
        
        self.stage_mtoward_button.clicked.connect(self.SM100_move_towardsample)
        self.stage_maway_button.clicked.connect(self.SM100_move_awaysample)
        
        self.stage_mmid_button.clicked.connect(self.SM100_move_mid)
        self.AskPositionButton.clicked.connect(self.AskPosition)
        """ parallel processing"""
        self.threadCount = QThreadPool.globalInstance().maxThreadCount()
        self.pool = QThreadPool.globalInstance()
        
        global stage_in_use
        stage_in_use = 0
        
    def initialise_stage(self):

        do_initialise = sub_initialise_stage(self.SMC100)
        self.pool.start(do_initialise)

    def disconnect_stage(self):
        self.SMC100.__del__()
        print('Stage disconnected. Power-off then -on for reconnection.')
    
    """ For up and down (camera view)
        Resolution of TRA-PP is 2.2 nm
        Travel range is 12 mm
        Minimum incremental motion is 100 nm
        self.RelativePos_box.value() is the multiple of 100 nm"""
        
    def SM100_move_left(self):
        axis = 1
        goto_pos = float(self.xystep_Box.value())*100/1e6
        do_move = sub_SM100_move(self.SMC100, str(axis), goto_pos)
        self.pool.start(do_move)
    
    def SM100_move_right(self):
        axis = 1
        goto_value = -float(self.xystep_Box.value())*100/1e6
        do_move = sub_SM100_move(self.SMC100, str(axis), goto_value)
        self.pool.start(do_move)
    
    def SM100_move_up(self):
        axis = 2
        goto_value = float(self.xystep_Box.value())*100/1e6
        do_move = sub_SM100_move(self.SMC100, str(axis), goto_value)
        self.pool.start(do_move)
    
    def SM100_move_down(self):
        axis = 2
        goto_value = -float(self.xystep_Box.value())*100/1e6
        do_move = sub_SM100_move(self.SMC100, str(axis), goto_value)
        self.pool.start(do_move)
        
    """ For z-movement (for focussing on sample)
        Resolution of LTA-HL is 7.4 nm
        Travel range is 25 mm
        Minimum incremental motion 50 nm
        self.RelativePos_box.value() is the multiple of 50 nm"""
    
    def SM100_move_towardsample(self):
        axis = 3
        goto_value = float(self.zstep_box.value())*50/1e6
        do_move = sub_SM100_move(self.SMC100, str(axis), goto_value)
        self.pool.start(do_move)
    
    def SM100_move_awaysample(self):
        axis = 3
        goto_value = -float(self.zstep_box.value())*50/1e6
        do_move = sub_SM100_move(self.SMC100, str(axis), goto_value)
        self.pool.start(do_move)
        
    def SM100_move_mid(self):
        do_move_mid = sub_SM100_move_mid(self.SMC100)
        self.pool.start(do_move_mid)
        
    def AskPosition(self):
        print(self.SMC100.get_position(1))
        print(self.SMC100.get_position(2))
        print(self.SMC100.get_position(3))
    
    def make_window(self):
        app = get_qt_app()
        self.show()
        app.exec_()
        return self
        
class sub_initialise_stage(QRunnable):
    def __init__(self, SMC100):
        super().__init__()
        self.SMC100 = SMC100

    def run(self):
        self.SMC100.reset_and_configure()
        self.SMC100.home()
        print('Stage initialised successful')
        print(self.SMC100.get_position(1))
        print(self.SMC100.get_position(2))
        print(self.SMC100.get_position(3))

class sub_SM100_move(QRunnable):
    def __init__(self, SMC100, axis, value):
        """ axis is a string and value is absolute position """
        super().__init__()
        self.SMC100 = SMC100
        self.axis = axis
        self.value = value
        global stage_in_use
        
    def run(self):
        global stage_in_use
        if stage_in_use > 0:
            print('Wait for stage release')
        while stage_in_use > 0:            
            time.sleep(0.5)
        print('Move')
        stage_in_use = 1
        pos_stage = self.SMC100.get_position(self.axis)
        goto_pos = float(pos_stage[0]) + self.value
        self.SMC100.move(goto_pos, self.axis, relative=False)
        pos_stage = self.SMC100.get_position(self.axis)
        print('Position axis ' + str(self.axis) + ': ' + str(pos_stage) + ' mm')
        stage_in_use = 0
        
class sub_SM100_move_mid(QRunnable):
    def __init__(self, SMC100):
        super().__init__()
        self.SMC100 = SMC100
        global stage_in_use

    def run(self):
        global stage_in_use
        if stage_in_use > 0:
            print('Wait for stage release')
        while stage_in_use > 0:            
            time.sleep(0.5)
        print('Move')
        stage_in_use = 1
        
        self.SMC100.move(6, '1', relative=False)
        self.SMC100.move(6, '2', relative=False)
        self.SMC100.move(6, '3', relative=False)
        print('Stage all axes moved to position 6 mm')
        stage_in_use = 0
        
""" For running GUI alone"""
#app = QtWidgets.QApplication(sys.argv)
#main = SMC100_window_object()
#main.show()
#sys.exit(app.exec_())








































