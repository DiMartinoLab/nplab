# -*- coding: utf-8 -*-
"""
Created on Jan 15 10:23:36 2019

@author: Hera
"""

import nplab

from nplab.instrument.spectrometer.seabreeze import OceanOpticsSpectrometer
from nplab.instrument.stage.SMC100 import SMC100
import nplab.utils.gui 
import nplab.datafile as datafile
from nplab.instrument.spectrometer import Spectrometer
from nplab.experiment import Experiment, ExperimentStopped
from nplab.utils.notified_property import DumbNotifiedProperty, NotifiedProperty
from nplab.utils.gui import QtCore, QtGui, uic, get_qt_app
from nplab.ui.ui_tools import UiTools
#import Rotation_Stage as RS
#import visa
import os
import time
import threading
import numpy as np
import pyqtgraph as pg
from nplab.ui.ui_tools import QtWidgets



#myOceanOptics = OceanOpticsSpectrometer(0)

class PL_mapping(Experiment, QtWidgets.QWidget, UiTools):
    # To use auto_connect_by_name name all widgets using _WidgetType, e.g. Vhigh_DoubleSpinBox
    # Then define DumbNotifiedProperty with the same name without _WidgetType, e.g. Vhigh
    
#    description = DumbNotifiedProperty("...")
#    single_Raman_spectrum_description = DumbNotifiedProperty('...')
    
    log_to_console = True
    live_Raman_spectrum_signal = QtCore.Signal(np.ndarray)
    live_darkfield_spectrum_signal = QtCore.Signal(np.ndarray)
    live_electronic_signal = QtCore.Signal(float, float, float)
    
    def __init__ (self, activeDatafile):
        super(PL_mapping, self).__init__()
        uic.loadUi('PL_mapping.ui', self)
        
###comment out software you are not going to use
#        self.initialise_smu() #Keithley, for electrical measurements
#        self.initialise_SmarAct_stage() #piezo stage for cantilever positioning##
        self.initialise_SMC100() #actuators for xy stage
#        self.initialise_OOSpectrometer() #for DF (white light) and PL (444nm laser)
#        self.initialise_shutter() #control box
#        self.initialise_Kandor() #Kymera, for Raman with 633nm or 785nm laser #jks68 19/10/2021
####end        
        self.radiantvoltages=None

        self.setup_plot_widgets()
        
        self.activeDatafile = activeDatafile

        self.auto_connect_by_name(self)
        
        self.live_darkfield_spectrum_signal.connect(self.update_darkfield_spectrum_plot)
        self.UploadFile.clicked.connect(self.processradiantfile)
        self.openstage.clicked.connect(self.open_SMC100_ui)
#        self.myArduino.shutterIN() #To ensure shutter is closed
        
#start of lab 3 experiment
    #TODO ALICE - go through this line by line
    def run(self, *args):
        # *args collects extra unnecessary arguments from qt
        self.voltages_data = []
        self.times_data = []
        self.currents_data = []
        self.capacitance_data = []
        stepwiseCounter = 1
        rampCounter = 1
        runningRampInterval = False
        try:
            activeDatagroup = self.create_data_group('scan_%d')
            activeDatagroup.attrs.create('description', str(self.description))
            if not self.mode_smuOnly.isChecked():
                if (self.mode_RamanOnly.isChecked() or self.mode_RamanAndDarkfield.isChecked() ):
                    self.RamanSpectrumImagePlotData = []
                    self.myKandor.AcquisitionMode = 1    # acquisition mode = 1 for single frame acquisition
                    self.myKandor.ReadMode = 3          # read mode = single track (reads centre_row +- num_rows). Output is one spectrum.
                    self.myKandor.set_camera_parameter('SingleTrack', self.centre_row, self.num_rows)
#                    self.myKandor.set_camera_parameter('Exposure', self.RamanIntegrationTime)
                    self.RamanWavelengths = self.myKandor.get_x_axis()
#                    self.RamanWavelengths = self.KandorControlUI.Andor.x_axis
                    self.RamanBackground = self.KandorControlUI.Andor.background
                    
                    activeDatagroup.attrs.create('Raman_Wavelengths', self.RamanWavelengths)
                    activeDatagroup.attrs.create('Raman_Background', self.RamanBackground)
                    activeDatagroup.attrs.create('Raman_Integration_Time', self.RamanIntegrationTime)
                    
                if (self.mode_DarkfieldOnly.isChecked() or self.mode_RamanAndDarkfield.isChecked() ):
                    self.darkfieldSpectrumImagePlotData = []
                    self.darkfieldWavelengths = self.OOspectrometer.read_wavelengths()
                    activeDatagroup.attrs.create('DarkfieldWavelengths', self.darkfieldWavelengths)
                    activeDatagroup.attrs.create('DarkfieldBackground', self.OOspectrometer.background)
                    activeDatagroup.attrs.create('DarkfieldReference', self.OOspectrometer.reference)
                    activeDatagroup.attrs.create('DarkfieldIntegrationTime', self.OOspectrometer.get_integration_time())
                    activeDatagroup.attrs.create('DarkfieldBackground_Reference', self.OOspectrometer.background_ref)
                    activeDatagroup.attrs.create('DarkfieldBackground_ReferenceIntegrationTime', self.OOspectrometer.background_int_ref)                    
                    #TEST BY THOMAS
                    #adding attributes for DF in-situ measurements
                    activeDatagroup.attrs.create('DarkfieldReferenceIntegrationTime', self.OOspectrometer.reference_int)
                    activeDatagroup.attrs.create('DarkfieldBackgroundIntegrationTime', self.OOspectrometer.background_int)
                if (self.mode_PLOnly.isChecked()):
                    self.darkfieldSpectrumImagePlotData = []
                    self.darkfieldWavelengths = self.OOspectrometer.read_wavelengths()
                    activeDatagroup.attrs.create('PLWavelengths', self.darkfieldWavelengths)
                    activeDatagroup.attrs.create('PLBackground', self.OOspectrometer.background)
                    activeDatagroup.attrs.create('PLIntegrationTime', self.OOspectrometer.get_integration_time())

            self.smu.src_voltage_range = float(self.Vrange_comboBox.currentText())
            self.smu.meas_current_range = float(self.Irange_comboBox.currentText())
            self.smu.src_voltage_limit = float(self.Vlimit_doubleSpinBox.value())
            self.smu.src_current_limit = float(self.Ilimit_doubleSpinBox.value())
            
            #TEST BY THOMAS
            #adding attributes for electrical measurements
            activeDatagroup.attrs.create('CurrentRange', self.smu.meas_current_range)
            activeDatagroup.attrs.create('CurrentLimit', self.smu.src_current_limit)
            try:
                activeDatagroup.attrs.create('VoltageProfilePts', self.radiantnopoints)
            except:
                print('No voltage profile uploaded!')
            
            activeVoltage = 0.0
            rampActiveVoltage = 0.0
            self.voltageRampSign = 1
            self.smu.output = 1
            self.smu.delay = 0
            t0 = time.time()
            print("----- Measurement started -----")
            n=0
            while True:
                self.smu.src_voltage = activeVoltage
                #spectrum = self.spectrometer.read_spectrum(bundle_metadata=True)
                if self.mode_smuOnly.isChecked():
                    self.acquireIVdatapoint(activeVoltage, t0, activeDatagroup)
                else:
                    
                    #TODO:
                    # for raman series: check mirror position
                    # put in DF configuration
                        # mirror in DF
                        # open shutter white light
                        # close andor shutter (not needed)
                        # take one DF spectrum
                    # put in raman configuration
                        # mirror in raman
                        # close shutter white light
                        # open andor shutter
                        # take raman, save
                    # put in DF configuration and save as DF final
                    
                    if (self.mode_RamanOnly.isChecked()):
                        RamanSpectrum_thread = threading.Thread(target = self.acquire_Raman_spectrum )    # acquiring Raman spectrum in new thread
                        RamanSpectrum_thread.start()
                        spectraActiveTime = (time.time()-t0)
                        activeDatagroup.append_dataset("RamanSpectra_times", spectraActiveTime)
                        while RamanSpectrum_thread.isAlive():       # collect current measurements while spectrum is being acquired
                            self.acquireIVdatapoint(activeVoltage, t0, activeDatagroup)
                        self.live_Raman_spectrum_signal.emit(self.RamanSpectrum)
                        activeDatagroup.append_dataset("RamanSpectrum", self.RamanSpectrum)
                    
                    if (self.mode_DarkfieldOnly.isChecked()):
                        DarkfieldSpectrum_thread = threading.Thread(target = self.acquire_darkfield_spectrum )    # acquiring Raman spectrum in new thread
                        DarkfieldSpectrum_thread.start()
                        spectraActiveTime = (time.time()-t0)
                        activeDatagroup.append_dataset("darkfieldSpectra_times", spectraActiveTime)
                        while DarkfieldSpectrum_thread.isAlive():       # collect current measurements while spectrum is being acquired
                            self.acquireIVdatapoint(activeVoltage, t0, activeDatagroup)
                        self.live_darkfield_spectrum_signal.emit(self.darkfieldSpectrum)
                        activeDatagroup.append_dataset("darkfieldSpectrum", self.darkfieldSpectrum)

                    if (self.mode_PLOnly.isChecked()):
                        DarkfieldSpectrum_thread = threading.Thread(target = self.acquire_darkfield_spectrum )    # acquiring Raman spectrum in new thread
                        DarkfieldSpectrum_thread.start()
                        spectraActiveTime = (time.time()-t0)
                        activeDatagroup.append_dataset("PLSpectra_times", spectraActiveTime)
                        while DarkfieldSpectrum_thread.isAlive():       # collect current measurements while spectrum is being acquired
                            self.acquireIVdatapoint(activeVoltage, t0, activeDatagroup)
                        self.live_darkfield_spectrum_signal.emit(self.darkfieldSpectrum)
                        activeDatagroup.append_dataset("PLSpectrum", self.darkfieldSpectrum)
                        
                    if (self.mode_RamanAndDarkfield.isChecked()):
                        # nothing for now
                        pass
                
                if self.hold.isChecked():       # if hold checkbox is checked do not change the voltage
                    pass
                else:                           # running in ramp or stepwise mode?
                    if self.RampON.isChecked():
                        if rampCounter < self.rampHoldNumber:
                            rampCounter += 1
                        else:                   # check if running a ramp with intermediate voltage value between steps
                            if (not self.RampIntermediateStep.isChecked()) or runningRampInterval:
                                if (rampActiveVoltage >= self.Vmax):    # if activeVoltage not in [Vmin,Vmax], e.g. after stepwise mode, reset it to Vmax or Vmin
                                    if (rampActiveVoltage > (self.Vmax + self.rampStep) ):
                                        rampActiveVoltage = self.Vmax + self.rampStep
                                    self.voltageRampSign = -1
                                elif (rampActiveVoltage <= self.Vmin):
                                    if (rampActiveVoltage < (self.Vmin - self.rampStep) ):
                                        rampActiveVoltage = self.Vmin - self.rampStep
                                    self.voltageRampSign = 1
                                rampActiveVoltage += self.voltageRampSign * self.rampStep
                                activeVoltage = rampActiveVoltage
                                runningRampInterval = False
                            elif self.RampIntermediateStep.isChecked() and (not runningRampInterval):
                                activeVoltage = self.rampIntermediateV
                                runningRampInterval = True
                            rampCounter = 1                            
                    elif self.RadiantFile.isChecked():
                        #self.smu_wait = self.radianttimedelay*0.001
                        time.sleep((1/(40*self.Vlow) - 1/40)) #DLCC
                        if type(self.radiantvoltages)==np.ndarray:
                            if n==self.radiantnopoints:
                                pass
                            activeVoltage = self.Vhigh * self.radiantvoltages[n]/100
                            n+=1
                            
                        else:
                            print('no file selected')
                            pass
                self.wait_or_stop(self.smu_wait)
                
        except ExperimentStopped:
            print("----- Measurement stopped -----")
        finally:
            self.activeDatafile.flush()
            self.smu.output = 0
            self.smu.src_voltage=0 #added by sunny to turn off voltage after experiement finished
            #self.AndorSpectrometer.light_shutter.open_shutter()

#end of lab 3 experiment


    def initialise_SMC100(self):
        self.SMC100=SMC100('COM1', (1,2,3))     
        
    def initialise_OOSpectrometer(self):
        self.OOspectrometer = OceanOpticsSpectrometer(0)
        

    def set_shamrock_wavelength(self):
        self.myKandor.kymera.SetWavelength(self.shamrockWavelength_spinBox.value())
        
    def setup_plot_widgets(self):
        self.electronics_plot = pg.PlotWidget()
        self.RamanSpectrum_plot = pg.PlotWidget()
        self.electronics_IVplot = pg.PlotWidget()
        self.RamanSpectrum_vs_time_plot = pg.PlotWidget()
        self.darkfieldSpectrum_plot = pg.PlotWidget()
        self.darkfieldSpectrum_vs_time_plot = pg.PlotWidget()
        self.capacitance_plot = pg.PlotWidget()
        self.replace_widget(self.plotGrid, self.plot1, self.electronics_plot)
        self.replace_widget(self.plotGrid, self.plot2, self.electronics_IVplot)
        self.replace_widget(self.plotGrid, self.plot3, self.RamanSpectrum_plot)
        self.replace_widget(self.plotGrid, self.plot4, self.RamanSpectrum_vs_time_plot)
#        self.replace_widget(self.plotGrid, self.plot5, self.darkfieldSpectrum_plot)
        self.replace_widget(self.plotGrid, self.plot5, self.capacitance_plot)
        self.replace_widget(self.plotGrid, self.plot6, self.darkfieldSpectrum_vs_time_plot)
#        # have to use pg.ImageItem() for an image plot. ImageItem is not a widget so it can't replace a PlotWidget directly,
#        # but it can be added as an ImageItem to a PlotWidget
        self.RamanImagePlotItem = pg.ImageItem()
        self.RamanSpectrum_vs_time_plot.addItem(self.RamanImagePlotItem)
        self.darkfieldImagePlotItem = pg.ImageItem()
        self.darkfieldSpectrum_vs_time_plot.addItem(self.darkfieldImagePlotItem)

       
    def update_darkfield_spectrum_plot(self, spectrum):
        self.darkfieldSpectrum_plot.plot(self.darkfieldWavelengths, spectrum, clear = True, pen = 'r')  # plot current darkfield spectrum in real time
        self.darkfieldSpectrumImagePlotData.append(np.nan_to_num(spectrum)[80:750])     # plot darkfield-time in range 400-900nm as image plot. ImageItem doesn't handle nan values, converting those to zero   
        self.darkfieldImagePlotItem.setImage(np.asarray(self.darkfieldSpectrumImagePlotData))
        

    def open_laser633_UI(self):
        self.laser633.show_gui()
        
    def open_SMC100_ui(self):
        self.SMC100.show_gui()
        
    def acquire_darkfield_spectrum(self):
        self.darkfieldSpectrum = np.asarray(self.OOspectrometer.read_processed_spectrum())
        
    def open_OO_spectrometer(self):
        self.gui_OOspectrometer= self.OOspectrometer.get_qt_ui()
        self.gui_OOspectrometer.show()
    
    
    def shutdown(self):
        self.activeDatafile.close()
        self.myKandor.cooler = False
        self.myKandor.close()
        self.OOspectrometer.shutdown_seabreeze()
        print('----Experiment ended----')
        #TODO - closes any UI windows
    
    def get_qt_ui(self):
        return self
    
if __name__ == '__main__':
    activeDatafile = nplab.current_datafile(working_directory = os.path.abspath(os.path.join(os.getcwd(),"../../..")))
    #working directory should now be C:\\Users\\<name>\\Documents    
    gui_activeDatafile = activeDatafile.get_qt_ui()
    gui_activeDatafile.show()
    
    experiment = PL_mapping(activeDatafile)
    experiment.show_gui()
    print('Done')
