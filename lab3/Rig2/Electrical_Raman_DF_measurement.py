
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

class smusetting():
    def __init__(self):
        self.Vmax = 0
        self.Vmin = 0
        self.rampStep = 0
        self.Vhigh = 0
        self.Vlow = 0
        self.smu_wait = 0
        self.rampStep = 0
        self.rampHoldNumber = 0
        self.rampIntermediateV = 0
        self.voltages_data = []
        self.times_data = []
        self.currents_data = []
        self.capacitance_data = []
        self.stepwiseCounter = 1
        self.rampCounter = 1
        self.runningRampInterval = False
    
class RadiantProfile():
    def __init__(self):
        self.radiantnopoints = 0 # first line describes how many points there are.
        self.radianttimedelay = 0 # second line describes time delay in ms between steps
        self.radiantvoltages = 0 # list of voltage steps
        
class Electrical_Raman_DF_measurement(Experiment, QtWidgets.QWidget, UiTools):
    #inherit QWidget methods
    def __init__(self, activeDatafile, Kymera_handle, Andor_handle,\
                 OOSpect_instance, smu_instance, myShutter, Ramanstate, threadCount, pool, \
                 ui_file = os.path.join(os.path.dirname(__file__),'Electrical_Raman_DF_measurement.ui'), \
                 pixel_number = 1024, pixel_width = 26, use_shifts = False, \
                     laser_wl = 782.5, white_shutter = None, settings_filepath = None, \
                     camera_index = None, parent = None):
        
        # make this window available for the main program
        super(Electrical_Raman_DF_measurement, self).__init__() 
        uic.loadUi(ui_file, self)
        
        self.kymera = Kymera_handle
        self.Andor = Andor_handle
        self.smu = smu_instance
        self.OOspectrometer = OOSpect_instance
        self.OOspectrometer.wl_in_nm = self.OOspectrometer.get_wavelengths()
        self.smusetting = smusetting()
        self.smu.smu_occupied = 0
        self.myShutter = myShutter
        self.Ramanstate = Ramanstate
        
        
        """ Prepare for file saving"""
        self.Datagroup = activeDatafile.create_group(name = 'SMU', auto_increment = True)
        self.Datagroup.attrs.create('Version Date', str('Date of latest edit of Raman Rig2 = 15-01-2024'))
                
        
        self.threadCount = threadCount
        self.pool = pool

        self.start_button.clicked.connect(self.startSMU)
        self.stop_button.clicked.connect(self.stopCapture)
        self.Measurementmode_button.clicked.connect(self.Switch_measurement_mode)
        self.smuSetLim_button.clicked.connect(self.smuSetLim)
        """Set default measurement mode to 'Raman only' """
        self.Measurementmode_button.setText('Raman only')
        self.smu.MeasurementMode = 'Raman only'
   
    def smuSetLim(self):
        self.smu.src_voltage_range = float(self.Vrange_comboBox.currentText())
        self.smu.meas_current_range = float(self.Irange_comboBox.currentText())
        self.smu.src_voltage_limit = float(self.Vlimit_doubleSpinBox.value())
        self.smu.src_current_limit = float(self.Ilimit_doubleSpinBox.value())
        
    def startSMU(self):
        
        # Set stopSMU = 0 as prevention
        self.smu.stopSMU = 0
        
        # check if smu is in use
        if self.smu.smu_occupied > 0:
            return
        
        # check if Andor is in use
        if self.Andor.Andor_occupied > 0:
            return
#        # check if OO is in use
#        if self.OOspectrometer.OO_occupied > 0:
#            return
        
        self.smusetting.Vmax = self.Vmax_doubleSpinBox.value()
        self.smusetting.Vmin = self.Vmin_doubleSpinBox.value()
        self.smusetting.rampStep = self.rampStep_doubleSpinBox.value()
        self.smusetting.Vhigh = self.Vhigh_doubleSpinBox.value()
        self.smusetting.Vlow = self.Vlow_doubleSpinBox.value()
        self.smusetting.smu_wait = self.smu_wait_doubleSpinBox.value()
        
        self.smusetting.rampStep = self.rampStep_doubleSpinBox.value()
        self.smusetting.rampHoldNumber = self.rampHoldNumber_spinbox.value()
        self.smusetting.rampIntermediateV = self.rampIntermediateV_doubleSpinBox.value()
    
    
        self.smusetting.voltages_data = []
        self.smusetting.times_data = []
        self.smusetting.currents_data = []
        self.smusetting.capacitance_data = []
        self.smusetting.stepwiseCounter = 1
        self.smusetting.rampCounter = 1
        self.smusetting.runningRampInterval = False
        self.smusetting.RampON = self.RampON
        self.smusetting.RadiantFile = self.RadiantFile
        self.smu.src_voltage_range = float(self.Vrange_comboBox.currentText())
        self.smu.meas_current_range = float(self.Irange_comboBox.currentText())
        self.smu.src_voltage_limit = float(self.Vlimit_doubleSpinBox.value())
        self.smu.src_current_limit = float(self.Ilimit_doubleSpinBox.value())
        
        
        if self.RadiantFile.isChecked():
            f = open("1000Hz_0.0010ramp 0.0010delay.txt", "r") #normal PUND with large gaps between pulses
            #f= open("constant-voltage-profile.txt", "r")
            #f = open("Thomas_shorter-PUND.txt", "r") #shorter PUND that reduces 0V gaps
            self.RadiantProfile.radiantnopoints = int(f.readline()[:-1]) #first line describes how many points there are.
            self.RadiantProfile.radianttimedelay = float(f.readline()[:-1]) #second line describes time delay in ms
            self.RadiantProfile.radiantvoltages = np.zeros(self.radiantnopoints)
            for x in range(self.RadiantProfile.radiantnopoints):
                self.RadiantProfile.radiantvoltages[x] = float(f.readline()[:-1])
            f.close()
            
        activeVoltage = 0.0
        rampActiveVoltage = 0.0
        self.voltageRampSign = 1
        
        self.smu.delay = 0
        t0 = time.time()
        print("----- Measurement started -----")
        n = 0
        
        #set starting voltage as 0 V
        self.smu.src_voltage = activeVoltage
        self.smu.output = 1
        
        self.Camerastate_output.setText('Capturing ...')
#                
        start_capture = sub_capture_and_save(self.Datagroup, \
                                             self.smu, \
                                             self.Andor, \
                                             self.kymera, \
                                             self.OOspectrometer, \
                                             self.Camerastate_output, \
                                             self.smusetting, \
                                             self.myShutter, \
                                             self.Ramanstate, \
                                             self.Filename_text, \
                                             self.Spectra_display, \
                                             self.Spectra_display_DF, \
                                             self.Voltage_display, \
                                             self.IV_display)
        self.smu.smu_occupied = 1
        self.pool.start(start_capture)

        
    def Switch_measurement_mode(self):
        available_modes = ['Raman only', 'Darkfield only', 'PL only', 'Raman/DF', 'SMU only']
        button_text = self.Measurementmode_button.text()
        if available_modes.index(button_text) == len(available_modes)-1:
            nextpos = 0
        else:
            nextpos = available_modes.index(button_text) + 1
        self.Measurementmode_button.setText(available_modes[nextpos])
        self.smu.MeasurementMode = available_modes[nextpos]
            
    def stopCapture(self):
        self.smu.stopSMU = 1
        self.smu.smu_occupied = 0
            
    def acquire_darkfield_spectrum(self):
#        self.darkfieldSpectrum = np.asarray(self.OOspectrometer.read_processed_spectrum())
        print('do acquire DF')
        self.darkfieldSpectrum = np.asarray(self.OOspectrometer.read_spectrum())
        
    def make_window(self):
        app = get_qt_app()
        self.show()
        app.exec_()
        return self


class sub_capture_and_save(QRunnable):
    def __init__(self, Datagroup, smu, Andor, kymera, OOspectrometer, \
                 Camerastate_output, smusetting, myShutter, Ramanstate, Filename_text, \
                 Spectra_display, Spectra_display_DF, Voltage_display, IV_display):
                 
        super().__init__()
        self.Datagroup = Datagroup
        self.smu = smu
        self.Andor = Andor
        self.kymera = kymera
        self.OOspectrometer = OOspectrometer
        self.Camerastate_output = Camerastate_output
        self.smusetting = smusetting
        self.Filename_text = Filename_text
        self.Spectra_display = Spectra_display
        self.Spectra_display_DF = Spectra_display_DF
        self.Voltage_display = Voltage_display
        self.IV_display = IV_display
        self.myShutter = myShutter
        self.Ramanstate = Ramanstate
        
    def run(self):

        """ Preparing to save txt file"""
        # take current time and date for txt save folder
        now = datetime.now()
        current_date = now.strftime("%y%m%d")
        # save each run into its own folder specified by time
        current_hour = now.strftime("%H%M")
        print('This run is saved under folder ' + str(current_hour))
        
        self.folder2save = os.path.join('C:\\Users\\Lab Di Martino\\Documents\\data\\' + self.OOspectrometer.CRSID)
        
        folders = []
        for root, dirs, files in os.walk(self.folder2save):
            for folder in dirs:
                if folder.startswith(current_date):
                    folders.append(folder)
        # create folder for the day if it is not made yet
        if not os.path.exists(self.folder2save + '\\' + str(folders[0])):
            os.makedirs(self.folder2save + '\\' + str(folders[0]))
        # create folder for the txt file save if it is not made yet    
        path2save = os.path.join(self.folder2save + '\\' + str(folders[0]) + '\\dataTxt')
        if not os.path.exists(path2save):
            os.makedirs(path2save)
        # create folder for this particular measurement
        path2save = os.path.join(self.folder2save + '\\' + str(folders[0]) + '\\dataTxt' + '\\' + current_hour)
        if not os.path.exists(path2save):
            os.makedirs(path2save)
        
        """ Create save name in h5-file"""
        if self.Filename_text.toPlainText() != 'File name ...':
            filename = self.Filename_text.toPlainText() + '_'
        else:
            filename = ''
        
        """ End preparing to save txt file"""
        
        """ Calculate Wavenumber for plotting only"""
        def WL2WN(wavelength):
                return 1/(self.Andor.laser_wl*1e-7) - 1/(wavelength*1e-7)
            
        def WN2WL(k):
            return 1e7/(1/(self.Andor.laser_wl*1e-7) - k)

        k = WL2WN(self.kymera.wavelength)
        """ END Calculate Wavenumber for plotting only"""
        
        t0 = time.time()
        self.Camerastate_output.setStyleSheet('background-color: rgb(255, 140, 140);')
        activeVoltage = self.smu.src_voltage
        voltageRampSign = 1
        
        
        
        """ clear the figures """
        self.Spectra_display.canvas.ax.cla()
        self.Spectra_display_DF.canvas.ax.cla()
        self.Voltage_display.canvas.ax.cla()
        self.IV_display.canvas.ax.cla()
        
        nscan = 0
        while self.smu.stopSMU < 1:
            
            if self.smusetting.RampON.isChecked():
                if self.smusetting.rampCounter < self.smusetting.rampHoldNumber:
                    self.smusetting.rampCounter += 1
                else:
                    if (activeVoltage + voltageRampSign * self.smusetting.rampStep >= self.smusetting.Vmax):
                        voltageRampSign = -1
                    elif (activeVoltage + voltageRampSign * self.smusetting.rampStep <= self.smusetting.Vmin):
                        voltageRampSign = 1
                    
                activeVoltage += voltageRampSign * self.smusetting.rampStep
                
                self.smu.src_voltage = activeVoltage
            
            spectraActiveTime = round(time.time()-t0, 2)
            
            print('doing measurement with ' + str(round(activeVoltage, 5)))
            
            """ Save current Voltage"""
            fullfile = path2save + '\\SMU_Voltage_' + filename + str(spectraActiveTime) + 's'
            outputfile = fullfile + '.txt'
            
            measuredCurrent = self.smu.read_current()
            
            afile = open(outputfile, 'w')
            saveData = np.array([[measuredCurrent], [activeVoltage]]).T
            np.savetxt(afile, saveData)
            afile.close()
            
            self.Voltage_display.canvas.ax.plot(nscan, activeVoltage, 'wo')
            self.Voltage_display.canvas.ax.set_xlabel('# scan', fontsize = 16, color = 'white')
            self.Voltage_display.canvas.ax.set_ylabel('Voltage (V)', fontsize = 16, color = 'white')
            self.Voltage_display.canvas.ax.tick_params(axis = "y", color = 'white',  direction = "in", right = False)
            self.Voltage_display.canvas.ax.tick_params(axis = "x", color = 'white', direction = "in", top = False)
            self.Voltage_display.canvas.fig.tight_layout()
            self.Voltage_display.canvas.draw()
            
            self.IV_display.canvas.ax.plot(activeVoltage, measuredCurrent*10e3, 'wo')
            self.IV_display.canvas.ax.set_xlabel('Voltage (V)', fontsize = 16, color = 'white')
            self.IV_display.canvas.ax.set_ylabel('Current (mA)', fontsize = 16, color = 'white')
            self.IV_display.canvas.ax.tick_params(axis = "y", color = 'white',  direction = "in", right = False)
            self.IV_display.canvas.ax.tick_params(axis = "x", color = 'white', direction = "in", top = False)
            self.IV_display.canvas.fig.tight_layout()
            self.IV_display.canvas.draw()
                
            if self.smu.MeasurementMode == 'Raman only' or self.smu.MeasurementMode == 'Raman/DF':
                
                # Close the DF shutter
                self.myShutter.query('B')
                # Flip mirror for Raman
                if int(self.Ramanstate.value()) == 1:
                    self.myShutter.query('F')
                elif int(self.Ramanstate.value()) == 2:
                    self.myShutter.query('E')
                # wait for flip mirror finished flipping
                time.sleep(3)
                
                self.Andor.Andor_occupied = 1
                imageArray, num_of_images, image_shape = self.Andor.capture()
                imageArray.reverse()
                spectrum = np.array(imageArray)
                
                # Release Andor occupation status
                self.Andor.Andor_occupied = 0
                
                """ Save Raman file"""
#                # Save as increment file name
#                i = 0
                
#                filenameR = filename + str(round(self.smu.src_voltage, 5)) + 'V_'
                
                # Save name with spectraActivetime
                fullfile = path2save + '\\Raman_SMU_' + filename + str(spectraActiveTime) + 's'
                outputfile = fullfile + '.txt'
            
            
                afile = open(outputfile, 'w')
                saveData = np.vstack(np.transpose((\
                                           np.append([0], self.kymera.wavelength), \
                                           np.append([round(self.Andor.Exposure, 2)], spectrum))))
                
                np.savetxt(afile, saveData)
                afile.close()
                """ Done save Raman file """
                
                
                self.Spectra_display.canvas.ax.pcolor([nscan, nscan+1], k, np.asarray([spectrum, spectrum*0]).T)
                self.Spectra_display.canvas.ax.set_xlabel('# scan', fontsize = 16, color = 'white')
                self.Spectra_display.canvas.ax.set_ylabel('Wavenumber (cm$^{-1}$)', fontsize = 16, color = 'white')
                self.Spectra_display.canvas.ax.tick_params(axis = "y", color = 'white',  direction = "in", right = False)
                self.Spectra_display.canvas.ax.tick_params(axis = "x", color = 'white', direction = "in", top = False)
                self.Spectra_display.canvas.fig.tight_layout()
                self.Spectra_display.canvas.draw()
            
            if self.smu.MeasurementMode == 'Darkfield only' or self.smu.MeasurementMode == 'Raman/DF': 
                if self.smu.MeasurementMode == 'Raman/DF':
                    while self.Andor.Andor_occupied == 1:
                        time.pause(1)
                    # Raman done
                    # Flip mirror for DF
                    if int(self.Ramanstate.value()) == 1:
                        self.myShutter.query('E')
                    elif int(self.Ramanstate.value()) == 2:
                        self.myShutter.query('F')
                    # wait for flip mirror finished flipping
                    time.sleep(3)
                
                # Open the DF shutter
                self.myShutter.query('A')
                time.sleep(3)
                self.OOspectrometer.OO_occupied = 1
                current_spec = self.OOspectrometer.read_spectrum()
                self.OOspectrometer.wl_in_nm
                # Release OO occupation status
                self.OOspectrometer.OO_occupied = 0
                
                """ Save DF file"""
                # Save as increment file name
                i = 0
                
                filenameR = filename + str(spectraActiveTime) + 'V_'
                fullfile = path2save + '\\DF_SMU_' + filenameR
                if os.path.exists(fullfile + '.txt'):
                    while os.path.exists(fullfile + '_' + str(i) + '.txt'):
                        i = i + 1
                    print('File DF nr. ' + str(i))
                    outputfile = fullfile + '_' + str(i) + '.txt'
                else:
                    outputfile = fullfile + '.txt'
            
                afile = open(outputfile, 'w')
                saveData = np.vstack(np.transpose((\
                                           np.append([0], self.OOspectrometer.wl_in_nm), \
                                           np.append([integrationtime], current_spec))))
                np.savetxt(afile, saveData)
                afile.close()
                """ Done save DF file """
                
                self.Spectra_display_DF.canvas.ax.pcolor([nscan, nscan+1], self.OOspectrometer.wl_in_nm, np.asarray([current_spec, current_spec*0]).T)
                self.Spectra_display_DF.canvas.ax.set_xlabel('# scan', fontsize = 16, color = 'white')
                self.Spectra_display_DF.canvas.ax.set_ylabel('Wavelength (nm)', fontsize = 16, color = 'white')
                self.Spectra_display_DF.canvas.ax.tick_params(axis = "y", color = 'white', direction = "in", right = False)
                self.Spectra_display_DF.canvas.ax.tick_params(axis = "x", color = 'white', direction = "in", top = False)
                self.Spectra_display_DF.canvas.fig.tight_layout()
                self.Spectra_display_DF.canvas.draw()

            if self.smu.MeasurementMode == 'PL only': 
                print('PL - no function yet')
                time.sleep(self.smusetting.smu_wait)
            
            
            if self.smu.MeasurementMode == 'SMU only': 
                time.sleep(self.smusetting.smu_wait)
                
            nscan += 1
            
        self.Camerastate_output.setText('Done!')
        self.Camerastate_output.setStyleSheet('background-color: rgb(180, 227, 255);')
        
        
        self.smu.output = 0
        self.smu.src_voltage = 0
        self.smu.smu_occupied = 0
        
















































        
        