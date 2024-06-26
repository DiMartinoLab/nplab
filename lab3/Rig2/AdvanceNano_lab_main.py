# -*- coding: utf-8 -*-
"""
Created on Jan 15 10:23:36 2019

@author: Hera
"""

import nplab
from nplab.instrument.camera.Andor.andor_sdk_rig2 import AndorBase
from nplab.instrument.spectrometer.Kymera import Kymera
from nplab.instrument.spectrometer.seabreeze import OceanOpticsSpectrometer, OceanOpticsControlUI
from nplab.instrument.electronics.keithley_2636b_smu import Keithley2636B as Keithley
# from nplab.instrument.stage.smaract_mcs import SmaractMCSSerial
from nplab.instrument.light_sources.matchbox_laser import MatchboxLaser
from nplab.instrument.stage.SMC100_lib_zstack import SMC100
import nplab.utils.gui 
import nplab.datafile as datafile
from nplab.instrument.spectrometer import Spectrometer
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
from PyQt5.QtCore import QTimer, QTime, Qt, QRunnable, QThreadPool

from nplab.ui.ui_tools import QtWidgets


import smaract.ctl as smaract_package
from lab3.Rig2.z_stack_window import z_stack_window_object
from lab3.Rig2.RamanSpectrometer import RamanSpectrometer
from lab3.Rig2.Electrical_Raman_DF_measurement import Electrical_Raman_DF_measurement
#from lab3.Rig2.MercuryController import MercuryController
from lab3.Rig2.SMC100StageControl_zstack import SMC100_window_object
from lab3.Rig2.Smaract_CryostatControl import Smaract_CryostatControl
from lab3.Rig2.xyz_stack_DF_Raman import xyz_stack_DF_Raman
#from lab3.OlympusCamera import OlympusCamera
from nplab.instrument.mercuryUSB.mercuryUSB import temperatureController as MiC_package
from lucam import Lucam
import sys
from nplab.instrument.serial_instrument import SerialInstrument

class AdvanceNano_lab_experiment(Experiment, QtWidgets.QWidget, UiTools):

    
    def __init__ (self, activeDatafile):
        super(AdvanceNano_lab_experiment, self).__init__()
        uic.loadUi('AdvanceNano_lab_main.ui', self)
        
        # which stage is enabled 'Cryostat stage' or 'Olympus stage' to load
        # into the z_stack UI
        self.stage_enabled = '1'
        
        """comment out software you are not going to use """
#        self.initialise_smu() #Keithley, for electrical measurements
#        self.initialise_SmarAct_stage() #piezo stage at cryostat       
#        self.initialise_MercuryControllers(truth_value = True) # Mercury controller iTC and iPS-M. do not initialise if truth_value is input as false       

        self.initialise_SMC100() #actuators for xy stage
        self.initialise_OOSpectrometer() #for DF (white light) and PL (444nm laser)
#        self.initialise_Kandor()
        
        self.initialise_shutter() #control box


#        self.initialise_camera() #Olympus camera
        

        self.radiantvoltages=None

        
        self.activeDatafile = activeDatafile
        self.singleRamanSpectraGroup = self.create_data_group('Single Raman spectra')
        self.z_stack_DataGroup = self.create_data_group('Data from z-stack window')

        self.auto_connect_by_name(self)
        
        self.openstage.clicked.connect(self.open_SMC100_ui)
        
        self.open_z_stack_window.clicked.connect(self.open_Dawn_z_stack_ui)

        self.OlympusCameraButton.clicked.connect(self.open_OlympusCamera)
        self.MercuryControllerButton.clicked.connect(self.open_MercuryController)
        self.FlipMirror_button.clicked.connect(self.flipthemirror)
        self.DFShutter_button.clicked.connect(self.DFshutter)
        self.Shutter_switch_button.clicked.connect(self.Raman_DF_switch)
        self.XYZStack_DF_Raman_button.clicked.connect(self.XYZStack_DF_Raman)
        self.open_Keithley_Andor_OO_UI_button.clicked.connect(self.open_electric_optic_measurement)
        # create general threadcount and pool for parallel measurements
        self.threadCount = QThreadPool.globalInstance().maxThreadCount()
        self.pool = QThreadPool.globalInstance()

    def initialise_Kandor(self):
        self.Kymera = Kymera()
        self.Andor = AndorBase()
        camera_index = None
        self.Andor.start(camera_index)
        print('Kymera and Andor initialised')     

    def initialise_smu(self):
        self.smu = Keithley.get_instance(address = 'USB0::0x05E6::0x2634::4454529::INSTR')
        self.smu.display = 0    # display current readings
#        self.smu.src_voltage_range = float(self.Vrange_comboBox.currentText())
#        self.smu.meas_current_range = float(self.Irange_comboBox.currentText())
#        self.smu.src_voltage_limit = float(self.Vlimit_doubleSpinBox.value())
#        self.smu.src_current_limit = float(self.Ilimit_doubleSpinBox.value())
        
    def initialise_SmarAct_stage(self):
        # Trung Nov. 2022
        # Code adapted from MSC2Example_DeviceInfo.py to control SmarAct stage
        # in rig2 for cryostat measurement
        try:
            buffer = smaract_package.FindDevices()
            if len(buffer) == 0:
                print("MCS2 no devices found.")
                sys.exit(1)
            locators = buffer.split("\n")
            for locator in locators:
                print("MCS2 available devices: {}".format(locator))
        except:
            print("MCS2 failed to find devices. Exit.")
            input()
            sys.exit(1)
        
        # Open the first MCS2 device from the list
        self.smaract_handle = smaract_package.Open(locators[0])
        # smaract_handle_enabled or SMC_handle_enabled to check
        # which stage z_stack_UI should load
        # check open_Dawn_z_stack_ui
        self.stage_enabled = 'Cryostat stage'
        print("MCS2 opened {}.".format(locators[0]))
        
    def initialise_MercuryControllers(self, truth_value):
        if truth_value ==True:
            self.MiTC_handle = MiC_package('COM8')
            self.MiPS_handle = MiC_package('COM10')
        elif truth_value == False:
            self.MiTC_handle, self.MiPS_handle = None, None
        
    def initialise_camera(self):
        self.camera_handle = Lucam()

    def initialise_SMC100(self):
        self.SMC100_handle = SMC100('COM1', (1,2,3))
        # smaract_handle_enabled or SMC_handle_enabled to check
        # which stage z_stack_UI should load
        # check open_Dawn_z_stack_ui
        self.stage_enabled = 'Olympus stage'
#    def initialise_z_stack(self):
#        self.z_stack_gui = z_stack()
        
    def initialise_shutter(self):
        self.myShutter = Arduino_shutter_class(port = 'COM9')
#        
#        self.myShutter.show_gui(blocking=False) # comment this out later. checking to make sure shutter works
        self.myShutter.termination_character = '\n'
        self.myShutter.port_settings = {
                    'baudrate':9600,
             #       'bytesize':serial.EIGHTBITS,
                    'timeout':2, #wait at most one second for a response
                    }
        self.myShutter.termination_character = '\r\n'
        
    def initialise_OOSpectrometer(self):
        self.OOspectrometer = OceanOpticsSpectrometer(0)


    
    # This is SmarAct stage for rotation
    def open_SmarAct_UI(self):
        delattr(self,'SmarAct_stage')
        self.SmarAct_stage = SmaractMCSSerial('COM6',3)
        self.SmarAct_stage.show_gui(blocking = False)

    def open_Andor_UI(self):
        if self.usercrsid_text.text() == 'CRSID of user ...':
            self.Andor.CRSID = 'xtn20'
        else:
            self.Andor.CRSID = 'self.usercrsid_text.text()'
            
##        self.KandorControlUI = self.myKandor.show_gui(block = False)
#        self.KandorControlUI = self.myKandor.get_control_widget()
#        self.KandorPreviewUI = self.myKandor.get_preview_widget()
#        self.KandorControlUI.show()
#        self.KandorPreviewUI.show()
        window_object_andor = RamanSpectrometer(activeDatafile = activeDatafile, \
                                          Kymera_handle = self.Kymera, \
                                          Andor_handle = self.Andor, \
                                          threadCount = self.threadCount, \
                                          pool = self.pool)
        self.this_ui_Andor = window_object_andor.make_window()
        self.open_Andor_UI_button.setEnabled(False)
        
    def open_electric_optic_measurement(self):
        
        if self.usercrsid_text.text() == 'CRSID of user ...':
            self.OOspectrometer.CRSID = 'xtn20'
        else:
            self.OOspectrometer.CRSID = 'self.usercrsid_text.text()'
            
        window_object_testandor = Electrical_Raman_DF_measurement(\
                                          activeDatafile = activeDatafile, \
                                          Kymera_handle = self.Kymera, \
                                          Andor_handle = self.Andor, \
                                          OOSpect_instance = self.OOspectrometer, \
                                          smu_instance = self.smu, \
                                          myShutter = self.myShutter, \
                                          Ramanstate = self.Ramanstate_box, \
                                          threadCount = self.threadCount, \
                                          pool = self.pool)
        self.this_ui_testandor = window_object_testandor.make_window()

    def open_laser633_UI(self):
        self.laser633.show_gui()
        
    def open_SMC100_ui(self):
        if self.stage_enabled == 'Cryostat stage':
            window_object_Smaract = Smaract_CryostatControl(Stage_instance = self.smaract_handle, \
                                          threadCount = self.threadCount, \
                                          pool = self.pool)
            self.this_ui_Smaract = window_object_Smaract.make_window()
        elif self.stage_enabled == 'Olympus stage':
            window_object_SMC = SMC100_window_object(SMC100_instance = self.SMC100_handle)
            self.this_ui_SMC = window_object_SMC.make_window()
            

    def open_Dawn_z_stack_ui(self):
        
        if self.usercrsid_text.text() == 'CRSID of user ...':
            self.OOspectrometer.CRSID = 'xtn20'
        else:
            self.OOspectrometer.CRSID = 'self.usercrsid_text.text()'
            
        print('THE Z-STACK BUTTON WAS CLICKED - '+self.OOspectrometer.get_model_name()+' was found.')
        # Create an instance of our class
        # if the temperature and magnet field controllers are not on
        #try:
         #   self.MiTC_handle
        #except:
         #   self.MiPS_handle = None
          #  self.MiTC_handle = None
            
        # check which stage is on to load 'Cryostat stage' or 'Olympus stage'
        print(self.stage_enabled)
        if self.stage_enabled == 'Cryostat stage':
            window_object_zstack = z_stack_window_object(activeDatafile=activeDatafile, \
                                                OOSpect_instance=self.OOspectrometer, \
                                                Stage_instance = self.smaract_handle, \
                                                MiTC_instance = self.MiTC_handle, \
                                                MiPS_instance = self.MiPS_handle, \
                                                stage_enabled = self.stage_enabled, \
                                                threadCount = self.threadCount, \
                                                pool = self.pool)
        elif self.stage_enabled == 'Olympus stage':
            window_object_zstack = z_stack_window_object(activeDatafile=activeDatafile, \
                                                OOSpect_instance=self.OOspectrometer, \
                                                Stage_instance = self.SMC100_handle, \
                                                MiTC_instance = None, MiPS_instance = None, \
                                                stage_enabled = self.stage_enabled, \
                                                threadCount = self.threadCount, \
                                                pool = self.pool)
        self.this_ui = window_object_zstack.make_window()
        del window_object_zstack
        self.open_z_stack_window.setEnabled(False)
    
    def XYZStack_DF_Raman(self):
        
        if self.usercrsid_text.text() == 'CRSID of user ...':
            self.OOspectrometer.CRSID = 'xtn20'
            self.Andor.CRSID = 'xtn20'
        else:
            self.OOspectrometer.CRSID = 'self.usercrsid_text.text()'
            self.Andor.CRSID = 'self.usercrsid_text.text()'
        
        # Create an instance of our class
        # if the temperature and magnet field controllers are not on
        #try:
         #   self.MiTC_handle
        #except:
         #   self.MiPS_handle = None
          #  self.MiTC_handle = None
            
        # check which stage is on to load 'Cryostat stage' or 'Olympus stage'
        print(self.stage_enabled)
        if self.stage_enabled == 'Cryostat stage':
            window_object_XYZStack_DF_Raman = xyz_stack_DF_Raman(activeDatafile=activeDatafile, \
                                                OOSpect_instance=self.OOspectrometer, \
                                                Stage_instance = self.smaract_handle, \
                                                MiTC_instance = self.MiTC_handle, \
                                                MiPS_instance = self.MiPS_handle, \
                                                stage_enabled = self.stage_enabled, \
                                                Kymera_handle = self.Kymera, \
                                                Andor_handle = self.Andor, \
                                                threadCount = self.threadCount, \
                                                pool = self.pool)
        elif self.stage_enabled == 'Olympus stage':
            window_object_XYZStack_DF_Raman = xyz_stack_DF_Raman(activeDatafile=activeDatafile, \
                                                OOSpect_instance=self.OOspectrometer, \
                                                Stage_instance = self.SMC100_handle, \
                                                MiTC_instance = None, MiPS_instance = None, \
                                                stage_enabled = self.stage_enabled, \
                                                Kymera_handle = self.Kymera, \
                                                Andor_handle = self.Andor, \
                                                threadCount = self.threadCount, \
                                                pool = self.pool)
        
        self.this_ui = window_object_XYZStack_DF_Raman.make_window()
        del window_object_XYZStack_DF_Raman
        self.XYZStack_DF_Raman_button.setEnabled(False)
        
    def open_OlympusCamera(self):
        print('Show me what you see camera!')
        window_object_camera = OlympusCamera(self.camera_handle) # Create an instance of our class
        self.this_ui = window_object_camera.make_window()
        
    def open_MercuryController(self):
        print('Run cool measurements!')
        window_object_MiTC = MercuryController(MiTC_handle = self.MiTC_handle) # Create an instance of our class
        self.this_ui = window_object_MiTC.make_window()

        
    def SmarAct_rotate_left(self):   
        self.SmarAct_stage.close()  
        self.SmarAct_stage = SmaractMCSSerial('COM6',4)
        SmarAct_rotation_steps = int( float(self.SmarAct_rotation_degrees_comboBox.currentText()) /0.0031 )   # from rough calibration 1 slip-stick step = 0.0031deg
        self.SmarAct_stage.slip_stick_move(3, steps = -SmarAct_rotation_steps, amplitude=1800, frequency=1000)   # channel 3 is rotation stage, default values for amplitude and frequency)
        
    def SmarAct_rotate_right(self):
        self.SmarAct_stage.close()
        self.SmarAct_stage = SmaractMCSSerial('COM6',4)
        SmarAct_rotation_steps = int( float(self.SmarAct_rotation_degrees_comboBox.currentText()) /0.0031 )
        self.SmarAct_stage.slip_stick_move(3, steps = SmarAct_rotation_steps, amplitude=1800, frequency=1000)

    def open_OO_spectrometer(self):
        self.gui_OOspectrometer= self.OOspectrometer.get_qt_ui()
        self.gui_OOspectrometer.show()

                
    def shutdown(self):
        self.activeDatafile.close()
        if hasattr(self, 'OOspectrometer'):
            self.OOspectrometer.shutdown_seabreeze()
            print('OO disconnected')
        if self.stage_enabled == 'Olympus stage':
            self.SMC100_handle.__del__()
            print('Olympus stage disconnected')
        # closes any UI windows
        win_list = QtWidgets.QApplication.allWindows()
        for w in win_list:
            w.close()
        print('----Experiment ended----')
    
    def get_qt_ui(self):
        return self
    
    def processradiantfile(self):
        #f= open("DLCC.txt", "r")
        f= open("1000Hz_0.0010ramp 0.0010delay.txt", "r") #normal PUND with large gaps between pulses
        #f= open("constant-voltage-profile.txt", "r")
        #f = open("Thomas_shorter-PUND.txt", "r") #shorter PUND that reduces 0V gaps
        self.radiantnopoints = int(f.readline()[:-1]) #first line describes how many points there are.
        self.radianttimedelay = float(f.readline()[:-1]) #second line describes time delay in ms
        self.radiantvoltages = np.zeros(self.radiantnopoints)
        for x in range(self.radiantnopoints):
            self.radiantvoltages[x] = float(f.readline()[:-1])
        f.close()
        
    
    def flipthemirror(self):
        if self.FlipMirror_button.text() == 'Flip mirror: 1':
            self.FlipMirror_button.setText('Flip mirror: 2')
            self.myShutter.query('E')
        elif self.FlipMirror_button.text() == 'Flip mirror: 2':
            self.FlipMirror_button.setText('Flip mirror: 1')
            self.myShutter.query('F')
        self.FlipMirror_button.setStyleSheet('background-color: rgb(0, 170, 127);')
    
    def DFshutter(self):
        if self.DFShutter_button.text() == 'DF shutter: Close':
            self.DFShutter_button.setText('DF shutter: Open')
            self.myShutter.query('A')
        elif self.DFShutter_button.text() == 'DF shutter: Open':
            self.DFShutter_button.setText('DF shutter: Close')
            self.myShutter.query('B')
        self.DFShutter_button.setStyleSheet('background-color: rgb(0, 170, 127);')
    
    def Raman_DF_switch(self):
        if self.Shutter_switch_button.text() == 'Switch to Raman':
            
#            print('to raman')
            # Close the DF shutter
            self.myShutter.query('B')
            self.DFShutter_button.setText('DF shutter: Close')
            # Flip mirror for Raman
            if int(self.Ramanstate_box.value()) == 1:
                self.FlipMirror_button.setText('Flip mirror: 1')
                self.myShutter.query('F')
            elif int(self.Ramanstate_box.value()) == 2:
                self.FlipMirror_button.setText('Flip mirror: 2')
                self.myShutter.query('E')
            self.Shutter_switch_button.setText('Switch to DF')
        
        elif self.Shutter_switch_button.text() == 'Switch to DF':
            
#            print('to df')
            # Open the DF shutter
            self.DFShutter_button.setText('DF shutter: Open')
            self.myShutter.query('A')
            # Flip mirror for DF
            if int(self.Ramanstate_box.value()) == 1:
                self.FlipMirror_button.setText('Flip mirror: 2')
                self.myShutter.query('E')
            elif int(self.Ramanstate_box.value()) == 2:
                self.FlipMirror_button.setText('Flip mirror: 1')
                self.myShutter.query('F')
            
            self.Shutter_switch_button.setText('Switch to Raman')  
        
        self.FlipMirror_button.setStyleSheet('background-color: rgb(0, 170, 127);')
        self.DFShutter_button.setStyleSheet('background-color: rgb(0, 170, 127);')
        self.Shutter_switch_button.setStyleSheet('background-color: rgb(0, 170, 127);')
        
class Arduino_shutter_class(SerialInstrument):
    '''Control class for tri shutter '''

    def __init__(self, port = None):
        '''Set up baudrate etc and recenters the stage to the center of it's range (50um)
        
        Args:
            port(int/str):  The port the device is connected 
                            to in any of the accepted serial formats
            
        '''
        self.termination_character = '\n'
        self.port_settings = {
                    'baudrate':9600,
             #       'bytesize':serial.EIGHTBITS,
                    'timeout':2, #wait at most one second for a response
                    }
        self.termination_character = '\r\n'
        SerialInstrument.__init__(self,port=port)

if __name__ == '__main__':
    activeDatafile = nplab.current_datafile(working_directory = os.path.abspath(os.path.join(os.getcwd(),"../../..")))
    #working directory should now be C:\\Users\\<name>\\Documents    
    gui_activeDatafile = activeDatafile.get_qt_ui()
    gui_activeDatafile.show()
    print(activeDatafile.keys())
    experiment = AdvanceNano_lab_experiment(activeDatafile)
    experiment.show_gui()

