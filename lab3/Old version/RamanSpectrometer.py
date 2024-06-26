
"""
Created on Wed Sep 13 16:35:28 2023

@author: Lab Di Martino
Trung - xtn20
September 2023
Interface and control for Andor spectrometer and camera

08-01-24 Script works
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

        
class RamanSpectrometer(Experiment, QtWidgets.QWidget, UiTools):
    #inherit QWidget methods
    def __init__(self, activeDatafile, Kymera_handle, Andor_handle,\
                 ui_file = os.path.join(os.path.dirname(__file__),'RamanSpectrometer.ui'), \
                 pixel_number = 1024, pixel_width = 26, use_shifts = False, \
                     laser_wl = 782.5, white_shutter = None, settings_filepath = None, \
                     camera_index = None, parent = None):
        
        # make this window available for the main program
        super(RamanSpectrometer, self).__init__() 
        uic.loadUi(ui_file, self)
        
        """ Initialise the spectrometer """
        self.kymera = Kymera_handle
        self.kymera.pixel_number = pixel_number
        self.kymera.pixel_width = pixel_width
        self.use_shifts = use_shifts
        self.laser_wl = laser_wl
        self.white_shutter = white_shutter
        self.pixel_number = pixel_number
        
        # set grating and central wavelength, uncomment while use with main GUI
        self.kymera.SetGrating(grating_num=1)
        self.kymera.SetWavelength(800)
#        self.metadata_property_names += ('slit_width', 'wavelengths')
        self.ImageFlip = 0
        global wavelength
        wavelength = self.get_x_axis()
        
        
        """ Initialise the camera """
        self.Andor = Andor_handle
        self.CurImage = None
        self.background = None
        self.backgrounded = False
        self.keep_shutter_open = False
        if settings_filepath is not None:
            self.load_params_from_file(settings_filepath)
        self.isAborted = False
        
        """ Buttons"""
        self.Centralwavelength_button.clicked.connect(self.set_centralwavelength)
        self.Exposure_button.clicked.connect(self.set_exposureTime)
        self.Capture_button.clicked.connect(self.CaptureSpectrum)
        self.Readmode_button.clicked.connect(self.Switch_read_mode)
        self.SetROI_button.clicked.connect(self.setROI)
        self.Abort_button.clicked.connect(self.Abort)
        self.Save_button.clicked.connect(self.Save)
        self.Repeat_button.clicked.connect(self.Capture_and_save)
        self.Live_button.clicked.connect(self.Live_mode)
        
        global Abort_live_mode
        Abort_live_mode = 0
        
        self.threadCount = QThreadPool.globalInstance().maxThreadCount()
        self.pool = QThreadPool.globalInstance()
        self.startGetTemp()
        # available_read_modes = ['FVB', 'Multi-track', 'Random track', 'Single track', 'Image']
        self.read_mode = 'Single track'
        self.acquisition_read_mode(self.read_mode)
        
        """ Prepare for file saving"""
        self.Datagroup = activeDatafile.create_group(name='AndorData', auto_increment=True)
        self.Datagroup.attrs.create('Version Date', str('Date of latest edit of Raman Rig2 = 15-09-2023'))
        
        """ update plot range"""
        self.update_figure_lim()
        self.SetLim_button.clicked.connect(self.update_figure_lim)
        
    # Spectrometer parts
    def get_x_axis(self):
        return self.kymera.GetCalibration()[::-1]
    
    def set_centralwavelength(self):
        self.kymera.SetWavelength(self.Centralwavelength_box.value())
        global wavelength
        wavelength = np.array(self.get_x_axis())
        self.update_figure_lim()
        
    # Camera parts
    def startGetTemp(self):
        print('Start logging')
        self.get_temperature()
        # Create a timer object
        timer = QTimer(self)
        # Add action to timer
        timer.timeout.connect(self.get_temperature)
        # Run the event each 30 seconds
        timer.start(30*1000)
        
    def get_temperature(self):
        read_temp = sub_get_temperature(self.Temperature_display, self.Andor)
        self.pool.start(read_temp)

    def acquisition_read_mode(self, readmode):
        available_modes = ['Single', 'Accumulate', 'Kinetic', 'Fast Kinetic']
        currentMode = 'Single'
        self.Andor.set_andor_parameter('AcquisitionMode', available_modes.index(currentMode) + 1)
        print('Acquisition mode set to Single')
        
        available_modes = ['FVB', 'Multi-track', 'Random track', 'Single track', 'Image']
        currentMode = readmode
        self.Andor.set_andor_parameter('ReadMode', available_modes.index(currentMode))
        print('Read mode set to ' + currentMode)
        
    def set_exposureTime(self):
        self.Andor.set_andor_parameter('Exposure', float(self.Exposure_box.value()))
        print('Exposure time set to ' + str(round(self.Andor.Exposure, 2)) + ' s')
        self.Exposure_box.setValue(round(self.Andor.Exposure, 2))
        
    def CaptureSpectrum(self):
        
        self.Camerastate_output.setText('Capturing ...')
        
        if self.read_mode == 'Single track':
            start_capture = sub_capture(self.Andor, self.Spectra_display, \
                                    self.Camerastate_output, self.laser_wl, \
                                    self.YLim1_box, self.YLim2_box)
        elif self.read_mode == 'Image':
            start_capture = sub_captureImage(self.Andor, self.Spectra_display, \
                                    self.Camerastate_output, self.laser_wl, \
                                    self.YLim1_box, self.YLim2_box)

        self.pool.start(start_capture)
        
    def Capture_and_save(self):
        self.Camerastate_output.setText('Capturing ...')
                
        start_capture = sub_capture_and_save(self.Datagroup, self.Andor, self.Spectra_display, \
                                             self.Camerastate_output, self.laser_wl, \
                                             self.Repeat_box.value(), self.Filename_text, \
                                             self.Description_text, self.YLim1_box, self.YLim2_box)
        self.pool.start(start_capture)
        
    def Abort(self):
       
        global Abort_live_mode
        Abort_live_mode = 1
        
    def Live_mode(self):
        self.Camerastate_output.setText('Capturing ...')
        self.Camerastate_output.setStyleSheet('background-color: rgb(255, 140, 140);')
        start_capture = sub_Live_mode(self.Andor, self.Spectra_display, \
                                             self.Camerastate_output, self.laser_wl,
                                             self.YLim1_box, self.YLim2_box)
        self.pool.start(start_capture)
        
        
    def Switch_read_mode(self):
        button_text = self.Readmode_button.text()
        if button_text == 'Image':
            self.read_mode = 'Single track'
            self.Readmode_button.setText('Single')
            self.acquisition_read_mode(self.read_mode)
        elif button_text == 'Single':
            self.read_mode = 'Image'
            self.Readmode_button.setText('Image')
            self.acquisition_read_mode(self.read_mode)
    
    def setROI(self):
        central_row = self.Centralrow_box.value()
        rows = self.Rows_box.value()
        print(self.read_mode)
        if self.read_mode == 'Image':
            self.Andor.set_andor_parameter('Image', 1, 1, 1, self.pixel_number, central_row - rows, central_row + rows)
        elif self.read_mode == 'Single track':
            self.Andor.set_andor_parameter('SingleTrack', central_row, rows*2)
    
    def Save(self):
        # for saving single measurement
        if self.Filename_text.toPlainText() != 'File name ...':
            filename = self.Filename_text.toPlainText()
        else:
            filename = 'Andor_data'
        
#        print(spectrum[0:3])
#        print(wavelength[0:3])
        activeRamandata = self.Datagroup.create_dataset(filename, data = spectrum)

        activeRamandata.attrs.create('Description', self.Description_text.toPlainText())        
        activeRamandata.attrs.create('Exposure_time', round(self.Andor.Exposure, 2))
        activeRamandata.attrs.create('wavelengths', wavelength)
        activeRamandata.attrs.create('background', 0)
        
        """ save data into text file"""
        now = datetime.now()
        current_date = now.strftime("%y%m%d")
        
        self.folder2save = os.path.join('C:\\Users\\Lab Di Martino\\Documents\\data\\xtn20')
        
        folders = []
        for root, dirs, files in os.walk(self.folder2save):
        
            for folder in dirs:
                if folder.startswith(current_date):
                    folders.append(folder)
        if not os.path.exists(self.folder2save + '\\' + str(folders[0])):
            os.makedirs(self.folder2save + '\\' + str(folders[0]))
            
        path2save = os.path.join(self.folder2save + '\\' + str(folders[0]) + '\\dataTxt')
        if not os.path.exists(path2save):
            os.makedirs(path2save)
        
        # Save as increment file name
        i = 0
        
        fullfile = path2save + '\\Raman_' + filename
        if os.path.exists(fullfile + '.txt'):
            while os.path.exists(fullfile + '_' + str(i) + '.txt'):
                print('file found')
                i = i + 1
            print(i)
            outputfile = fullfile + '_' + str(i) + '.txt'
        else:
            outputfile = fullfile + '.txt'
            
        afile = open(outputfile, 'w')
        saveData = np.vstack(np.transpose((wavelength, spectrum)))
        np.savetxt(afile, saveData)
        afile.close()


        """ END: save data into text file"""
       
    def update_figure_lim(self):
        
        self.laser_wl = float(self.LaserWavelength_value.value())
        self.Spectra_display.fontsize = 16
        def WL2WN(wavelength):
            return 1/(self.laser_wl*1e-7) - 1/(wavelength*1e-7)
        
        def WN2WL(k):
            return 1e7/(1/(self.laser_wl*1e-7) - k)
        
        self.Spectra_display.canvas.ax.set_ylabel('Intensity (counts)', fontsize = self.Spectra_display.fontsize, color = 'white')
        self.Spectra_display.canvas.ax.set_xlabel('Wavelength (nm)', fontsize = self.Spectra_display.fontsize, color = 'white')
        self.Spectra_display.canvas.ax.tick_params(axis = "y", direction = "in")
        self.Spectra_display.canvas.ax.tick_params(axis = "x", direction = "in", top = False)
        self.Spectra_display.canvas.ax.set_xlim([min(wavelength), max(wavelength)])
        
        self.Spectra_display.minY = min(float(self.YLim1_box.value()), float(self.YLim2_box.value()))
        self.Spectra_display.maxY = max(float(self.YLim1_box.value()), float(self.YLim2_box.value()))
        if self.Spectra_display.minY < self.Spectra_display.maxY:
            self.Spectra_display.canvas.ax.set_ylim([self.Spectra_display.minY, self.Spectra_display.maxY])
        
        """ Display WN axis"""
        
        self.Spectra_display.minX = min(float(self.XLim1_box.value()), float(self.XLim2_box.value()))
        self.Spectra_display.maxX = max(float(self.XLim1_box.value()), float(self.XLim2_box.value()))
        if self.Spectra_display.minX == self.Spectra_display.maxX:
            self.Spectra_display.minX = min(wavelength)
            self.Spectra_display.maxX = max(wavelength)
            
        self.Spectra_display.canvas.ax.set_xlim([self.Spectra_display.minX, self.Spectra_display.maxX])
        
        self.Spectra_display.canvas.ax2.set_xlabel('Wavenumber (cm$^{-1}$)', fontsize = self.Spectra_display.fontsize, color = 'white')
#        self.Spectra_display.canvas.ax2.set_xlim([self.Spectra_display.minX, self.Spectra_display.maxX])
        
        kmin = WL2WN(self.Spectra_display.minX)
        kmax = WL2WN(self.Spectra_display.maxX)
        kticksn = list(range(0, round(kmin), -200))
        kticksn.reverse()
        kticksp = list(range(0, round(kmax), 200))
        kticksp1 = list(range(0, 250, 50))
        kticksp2 = list(range(400, round(kmax), 200))
        kticksp = [*kticksp1, *kticksp2]
        kticks = [*kticksn, *kticksp]
        kticksInNm = WN2WL(np.array(kticks))
        kticksstr = []
        for i in kticks:
            kticksstr.append(str(i))
        self.Spectra_display.canvas.ax2.set_xticks(kticksInNm)
        self.Spectra_display.canvas.ax2.set_xticklabels(kticksstr)
        self.Spectra_display.canvas.ax2.set_xlim([self.Spectra_display.minX, self.Spectra_display.maxX])
        """ Display WN axis"""
        
        for label in (self.Spectra_display.canvas.ax.get_xticklabels() +  \
                      self.Spectra_display.canvas.ax.get_yticklabels() + \
                      self.Spectra_display.canvas.ax2.get_xticklabels()):
            label.set_fontsize(self.Spectra_display.fontsize)
        self.Spectra_display.canvas.fig.tight_layout()
        self.Spectra_display.canvas.draw()
        
        
    def make_window(self):
        app = get_qt_app()
        self.show()
        app.exec_()
        return self

#app = QtWidgets.QApplication(sys.argv)
#main = RamanSpectrometer()
#main.show()
#sys.exit(app.exec_())


class sub_get_temperature(QRunnable):
    def __init__(self, Temperature_display, Andor):
        super().__init__()
        self.Temperature_display = Temperature_display
        self.Andor = Andor
        
    def run(self):
        self.Temperature_display.display(float(self.Andor.CurrentTemperature))
        
class sub_capture(QRunnable):
    def __init__(self, Andor, Spectra_display, Camerastate_output, laser_wl, YLim1_box, YLim2_box):
        super().__init__()
        self.Andor = Andor
        self.Spectra_display = Spectra_display
        self.Camerastate_output = Camerastate_output
        self.laser_wl = laser_wl
        self.YLim1_box = YLim1_box
        self.YLim2_box = YLim2_box
        
    def run(self):
        global spectrum
        global wavelength
        
        self.Camerastate_output.setStyleSheet('background-color: rgb(255, 140, 140);')
        
        imageArray, num_of_images, image_shape = self.Andor.capture()
        imageArray.reverse()
        
        spectrum = np.array(imageArray)

        
        self.Spectra_display.canvas.ax.cla()

        def WL2WN(wavelength):
            return 1/(self.laser_wl*1e-7) - 1/(wavelength*1e-7)
        
        def WN2WL(k):
            return 1e7/(1/(self.laser_wl*1e-7) - k)

        p = self.Spectra_display.canvas.ax.plot(wavelength, spectrum)
        self.Spectra_display.canvas.ax.set_ylabel('Intensity (counts)', fontsize = 16, color = 'white')
        self.Spectra_display.canvas.ax.set_xlabel('Wavelength (nm)', fontsize = 16, color = 'white')
        self.Spectra_display.canvas.ax.tick_params(axis = "y", direction = "in")
        self.Spectra_display.canvas.ax.tick_params(axis = "x", direction = "in", top = False)
            
        self.Spectra_display.minY = min(float(self.YLim1_box.value()), float(self.YLim2_box.value()))
        self.Spectra_display.maxY = max(float(self.YLim1_box.value()), float(self.YLim2_box.value()))
        
        if self.Spectra_display.minY < self.Spectra_display.maxY:
            self.Spectra_display.canvas.ax.set_ylim([self.Spectra_display.minY, self.Spectra_display.maxY])
        self.Spectra_display.canvas.ax.set_xlim([self.Spectra_display.minX, self.Spectra_display.maxX])
        self.Spectra_display.canvas.ax2.set_xlim([self.Spectra_display.minX, self.Spectra_display.maxX])

#        
        for label in (self.Spectra_display.canvas.ax.get_xticklabels() +  \
                      self.Spectra_display.canvas.ax.get_yticklabels() + \
                      self.Spectra_display.canvas.ax2.get_xticklabels()):
            label.set_fontsize(self.Spectra_display.fontsize)
        self.Spectra_display.canvas.fig.tight_layout()
        self.Spectra_display.canvas.draw()
        
        self.Camerastate_output.setText('Done!')
        self.Camerastate_output.setStyleSheet('background-color: rgb(180, 227, 255);')
        
        """ unused
        #        self.Spectra_display.canvas.ax.set_ylabel('Intensity (counts)', fontsize = 16)
#        self.Spectra_display.canvas.ax.set_xlabel('Wavelength (nm)', fontsize = 16)
#        self.Spectra_display.canvas.ax.tick_params(axis = "y", direction = "in")
#        self.Spectra_display.canvas.ax.tick_params(axis = "x", direction = "in", top = False)
#        self.Spectra_display.canvas.ax.set_xlim([min(wavelength), max(wavelength)])
#        
#        minY = min(float(self.YLim1_box.value()), float(self.YLim2_box.value()))
#        maxY = max(float(self.YLim1_box.value()), float(self.YLim2_box.value()))
#        self.Spectra_display.canvas.ax.set_ylim([minY, maxY])
#        
#        self.Spectra_display.canvas.ax2.set_xlabel('Wavenumber (cm$^{-1}$)', fontsize = 16)
#        self.Spectra_display.canvas.ax2.set_xlim(\
#                        [min(wavelength), max(wavelength)])
#        
#        kmin = WL2WN(min(wavelength))
#        kmax = WL2WN(max(wavelength))
#        kticksn = list(range(0, round(kmin), -200))
#        kticksn.reverse()
#        kticksp = list(range(0, round(kmax), 200))
#        kticks = [*kticksn, *kticksp]
#        kticksInNm = WN2WL(np.array(kticks))
#        kticksstr = []
#        for i in kticks:
#            kticksstr.append(str(i))
#        self.Spectra_display.canvas.ax2.set_xticks(kticksInNm)
#        self.Spectra_display.canvas.ax2.set_xticklabels(kticksstr)

        """
        
        
class sub_capture_and_save(QRunnable):
    def __init__(self, Datagroup, Andor, Spectra_display, Camerastate_output, \
                 laser_wl, Repeat_box_value, Filename_text, Description_text, YLim1_box, YLim2_box):
        
        super().__init__()
        self.Datagroup = Datagroup
        self.Andor = Andor
        self.Spectra_display = Spectra_display
        self.Camerastate_output = Camerastate_output
        self.laser_wl = laser_wl
        self.Repeat_box_value = Repeat_box_value
        self.Filename_text = Filename_text
        self.Description_text = Description_text
        self.YLim1_box = YLim1_box
        self.YLim2_box = YLim2_box
        
    def run(self):
        global spectrum
        global wavelength
        
        i = 0
        while i <= self.Repeat_box_value:
            print(i)
            self.Camerastate_output.setText('Capturing ... ' + str(i))
            self.Camerastate_output.setStyleSheet('background-color: rgb(255, 140, 140);')
            imageArray, num_of_images, image_shape = self.Andor.capture()
            imageArray.reverse()
            
            spectrum = np.array(imageArray)
           
            self.Spectra_display.canvas.ax.cla()

            
            wavelength2 = np.array(wavelength)
            def WL2WN(wavelength):
                return 1/(self.laser_wl*1e-7) - 1/(wavelength*1e-7)
            
            def WN2WL(k):
                return 1e7/(1/self.laser_wl*1e-7 - k)
            
            def WN2WL(k):
                return 1e7/(1/(self.laser_wl*1e-7) - k)

            k = WL2WN(wavelength2)
    
            p = self.Spectra_display.canvas.ax.plot(wavelength, spectrum)
            self.Spectra_display.canvas.ax.set_ylabel('Intensity (counts)', fontsize = 16, color = 'white')
            self.Spectra_display.canvas.ax.set_xlabel('Wavelength (nm)', fontsize = 16, color = 'white')
            self.Spectra_display.canvas.ax.tick_params(axis = "y", direction = "in")
            self.Spectra_display.canvas.ax.tick_params(axis = "x", direction = "in", top = False)
            
            minY = min(float(self.YLim1_box.value()), float(self.YLim2_box.value()))
            maxY = max(float(self.YLim1_box.value()), float(self.YLim2_box.value()))
            if self.Spectra_display.minY < self.Spectra_display.maxY:
                self.Spectra_display.canvas.ax.set_ylim([self.Spectra_display.minY, self.Spectra_display.maxY])
            self.Spectra_display.canvas.ax.set_xlim([self.Spectra_display.minX, self.Spectra_display.maxX])
            self.Spectra_display.canvas.ax2.set_xlim([self.Spectra_display.minX, self.Spectra_display.maxX])
    
            for label in (self.Spectra_display.canvas.ax.get_xticklabels() +  \
                          self.Spectra_display.canvas.ax.get_yticklabels() + \
                          self.Spectra_display.canvas.ax2.get_xticklabels()):
                label.set_fontsize(self.Spectra_display.fontsize)
            self.Spectra_display.canvas.fig.tight_layout()
            self.Spectra_display.canvas.draw()
            
    
            """ save file """
            if self.Filename_text.toPlainText() != 'File name ...':
                filename = self.Filename_text.toPlainText()
            else:
                filename = 'Andor_data'
                
            activeRamandata = self.Datagroup.create_dataset(filename, data = spectrum)
    
            activeRamandata.attrs.create("Description", self.Description_text.toPlainText())        
            activeRamandata.attrs.create("Exposure_time", round(self.Andor.Exposure, 2))
            activeRamandata.attrs.create("wavelengths", wavelength)
            activeRamandata.attrs.create("background", 0)
            i += 1
            
        self.Camerastate_output.setText('Done!')
        self.Camerastate_output.setStyleSheet('background-color: rgb(180, 227, 255);')





class sub_Live_mode(QRunnable):
    def __init__(self, Andor, Spectra_display, Camerastate_output, \
                 laser_wl, YLim1_box, YLim2_box):
        
        super().__init__()

        self.Andor = Andor
        self.Spectra_display = Spectra_display
        self.Camerastate_output = Camerastate_output
        self.laser_wl = laser_wl
        self.YLim1_box = YLim1_box
        self.YLim2_box = YLim2_box
        
    def run(self):
        global spectrum
        global wavelength
        global Abort_live_mode
        

        while Abort_live_mode < 1:

            self.Camerastate_output.setText('Capturing ... ')
            self.Camerastate_output.setStyleSheet('background-color: rgb(255, 140, 140);')
            imageArray, num_of_images, image_shape = self.Andor.capture()
            imageArray.reverse()
            
            spectrum = np.array(imageArray)
           
            self.Spectra_display.canvas.ax.cla()
            
            wavelength2 = np.array(wavelength)
            def WL2WN(wavelength):
                return 1/(self.laser_wl*1e-7) - 1/(wavelength*1e-7)
            
            def WN2WL(k):
                return 1e7/(1/(self.laser_wl*1e-7) - k)

            k = WL2WN(wavelength2)
    
            p = self.Spectra_display.canvas.ax.plot(wavelength, spectrum)
            self.Spectra_display.canvas.ax.set_ylabel('Intensity (counts)', fontsize = 16, color = 'white')
            self.Spectra_display.canvas.ax.set_xlabel('Wavelength (nm)', fontsize = 16, color = 'white')
            self.Spectra_display.canvas.ax.tick_params(axis = "y", direction = "in")
            self.Spectra_display.canvas.ax.tick_params(axis = "x", direction = "in", top = False)
            
            minY = min(float(self.YLim1_box.value()), float(self.YLim2_box.value()))
            maxY = max(float(self.YLim1_box.value()), float(self.YLim2_box.value()))
            if self.Spectra_display.minY < self.Spectra_display.maxY:
                self.Spectra_display.canvas.ax.set_ylim([self.Spectra_display.minY, self.Spectra_display.maxY])
            self.Spectra_display.canvas.ax.set_xlim([self.Spectra_display.minX, self.Spectra_display.maxX])
            self.Spectra_display.canvas.ax2.set_xlim([self.Spectra_display.minX, self.Spectra_display.maxX])
    
            for label in (self.Spectra_display.canvas.ax.get_xticklabels() +  \
                          self.Spectra_display.canvas.ax.get_yticklabels() + \
                          self.Spectra_display.canvas.ax2.get_xticklabels()):
                label.set_fontsize(self.Spectra_display.fontsize)
            self.Spectra_display.canvas.fig.tight_layout()
            self.Spectra_display.canvas.draw()
            
    
        Abort_live_mode = 0
        self.Camerastate_output.setText('Done!')
        self.Camerastate_output.setStyleSheet('background-color: rgb(180, 227, 255);')
        

class sub_captureImage(QRunnable):
    def __init__(self, Andor, Spectra_display, Camerastate_output, laser_wl, YLim1_box, YLim2_box):
        super().__init__()
        self.Andor = Andor
        self.Spectra_display = Spectra_display
        self.Camerastate_output = Camerastate_output
        self.laser_wl = laser_wl
        self.YLim1_box = YLim1_box
        self.YLim2_box = YLim2_box
        
    def run(self):
        global spectrum
        global wavelength
        
        self.Camerastate_output.setStyleSheet('background-color: rgb(255, 140, 140);')
        
        imageArray, num_of_images, image_shape = self.Andor.capture()
        imageArray.reverse()
        
        """ imageArray is 1D array with length = x_pixel*y_pixel of the camera.
        We need to decompose the imageArray into 2D array"""
        
        image2D = np.array([ [0] * image_shape[1] for _ in range(int(image_shape[0]))])

        for i in np.arange(image_shape[0]):
            for j in np.arange(image_shape[1]):
                image2D[i][j] = imageArray[i * image_shape[1] + j]
        
        """ End of decomposing"""
        
        # Create x and y axis as number of pixel
        
        real_y = np.arange(image_shape[0])
        real_x = wavelength
        
        dx = (real_x[1]-real_x[0])/2.
        dy = (real_y[1]-real_y[0])/2.
        extent = [real_x[0]-dx, real_x[-1]+dx, real_y[0]-dy, real_y[-1]+dy]
        
        self.Spectra_display.canvas.ax.cla()
        self.Spectra_display.canvas.ax.imshow(image2D, extent = extent, aspect='auto')
        self.Spectra_display.canvas.ax.set_ylabel('Pixel (nr.)', fontsize = 16, color = 'white')
        self.Spectra_display.canvas.ax.set_xlabel('Wavelength (nm)', fontsize = 16, color = 'white')
        self.Spectra_display.canvas.ax.tick_params(axis = "y", direction = "in")
        self.Spectra_display.canvas.ax.tick_params(axis = "x", direction = "in", top = False)
            
        self.Spectra_display.canvas.ax.set_ylim([0, image_shape[0]])
        
        self.Spectra_display.canvas.ax.set_xlim([self.Spectra_display.minX, self.Spectra_display.maxX])
        self.Spectra_display.canvas.ax2.set_xlim([self.Spectra_display.minX, self.Spectra_display.maxX])

#        
        for label in (self.Spectra_display.canvas.ax.get_xticklabels() +  \
                      self.Spectra_display.canvas.ax.get_yticklabels() + \
                      self.Spectra_display.canvas.ax2.get_xticklabels()):
            label.set_fontsize(self.Spectra_display.fontsize)
        self.Spectra_display.canvas.fig.tight_layout()
        self.Spectra_display.canvas.draw()
        
        self.Camerastate_output.setText('Done!')
        self.Camerastate_output.setStyleSheet('background-color: rgb(180, 227, 255);')











