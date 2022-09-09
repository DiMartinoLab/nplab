# -*- coding: utf-8 -*-
"""
Created on Jan 15 10:23:36 2019

@author: Hera
"""

import nplab

from nplab.instrument.spectrometer.kandor import Kandor
from nplab.instrument.spectrometer.seabreeze import OceanOpticsSpectrometer
from nplab.instrument.electronics.keithley_2636b_smu import Keithley2636B as Keithley
#from nplab.instrument.stage.smaract_mcs import SmaractMCSSerial
from nplab.instrument.shutter.Arduino_ttl_shutter import Arduino_tri_shutter as shutter
from nplab.instrument.light_sources.matchbox_laser import MatchboxLaser

from nplab.instrument.stage.SMC100 import SMC100
from lab3.Dawn_z_stack import z_stack
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

class Lab3_experiment(Experiment, QtWidgets.QWidget, UiTools):
    # To use auto_connect_by_name name all widgets using _WidgetType, e.g. Vhigh_DoubleSpinBox
    # Then define DumbNotifiedProperty with the same name without _WidgetType, e.g. Vhigh
    Vmax = DumbNotifiedProperty(1.0)
    Vmin = DumbNotifiedProperty(-1.0)
    rampStep = DumbNotifiedProperty(0.1)
    Vhigh = DumbNotifiedProperty(1.0)
    Vlow = DumbNotifiedProperty(1.0)
    smu_wait = DumbNotifiedProperty(0.0)
    RamanIntegrationTime = DumbNotifiedProperty(1.0)
    laser633_power = DumbNotifiedProperty(0)
    centre_row = DumbNotifiedProperty(128)
    num_rows = DumbNotifiedProperty(128)
    stepwiseHoldNumber = DumbNotifiedProperty(3)
    rampHoldNumber = DumbNotifiedProperty(1)
    rampIntermediateV = DumbNotifiedProperty(0.0)
    
    description = DumbNotifiedProperty("description")
    single_Raman_spectrum_description = DumbNotifiedProperty('single spectrum description')
    
    log_to_console = True
    live_Raman_spectrum_signal = QtCore.Signal(np.ndarray)
    live_darkfield_spectrum_signal = QtCore.Signal(np.ndarray)
    live_electronic_signal = QtCore.Signal(float, float, float)
    
    def __init__ (self, activeDatafile):
        super(Lab3_experiment, self).__init__()
        uic.loadUi('lab3_interface_kymera.ui', self)
        
###comment out software you are not going to use
#        self.initialise_smu() #Keithley, for electrical measurements
#        self.initialise_SmarAct_stage() #piezo stage for cantilever positioning##
        self.initialise_SMC100() #actuators for xy stage
        self.initialise_OOSpectrometer() #for DF (white light) and PL (444nm laser)
#        self.initialise_shutter() #control box
#        self.initialise_Kandor() #Kymera, for Raman with 633nm or 785nm laser #jks68 19/10/2021
####end        
        self.radiantvoltages=None

        self.setup_plot_widgets()
        
        self.activeDatafile = activeDatafile
        self.singleRamanSpectraGroup = self.create_data_group('Single Raman spectra')

        self.auto_connect_by_name(self)
        
        self.live_Raman_spectrum_signal.connect(self.update_Raman_spectrum_plot)
        self.live_electronic_signal.connect(self.update_electronic_plot)
        self.live_darkfield_spectrum_signal.connect(self.update_darkfield_spectrum_plot)
#        self.andor_cooler_checkBox.toggled.connect(self.andor_cooler)
        self.UploadFile.clicked.connect(self.processradiantfile)
        self.openstage.clicked.connect(self.open_SMC100_ui)
        self.z_stack.clicked.connect(self.open_Dawn_z_stack_ui)
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

    def initialise_Kandor(self):
        self.myKandor = Kandor()
        print('Kymera initialised')     

    def initialise_smu(self):
        self.smu = Keithley.get_instance(address = 'USB0::0x05E6::0x2634::4454529::INSTR')
        self.smu.display = 0    # display current readings
        self.smu.src_voltage_range = float(self.Vrange_comboBox.currentText())
        self.smu.meas_current_range = float(self.Irange_comboBox.currentText())
        self.smu.src_voltage_limit = float(self.Vlimit_doubleSpinBox.value())
        self.smu.src_current_limit = float(self.Ilimit_doubleSpinBox.value())
        
    def initialise_SmarAct_stage(self):
        self.SmarAct_stage = SmaractMCSSerial('COM6',3)

    def initialise_SMC100(self):
        self.SMC100=SMC100('COM1', (1,2,3))     
        
    def initialise_shutter(self):
        self.myShutter = shutter(port = 'COM4')
        time.sleep(5)
        self.myShutter.show_gui(blocking=False) # comment this out later. checking to make sure shutter works
        
    def initialise_OOSpectrometer(self):
        self.OOspectrometer = OceanOpticsSpectrometer(0)
        
    def initialise_laser633(self):
        self.laser633 = MatchboxLaser("COM3")
    
#jks68 20/10/2021 start #CHECK WHETHER IT WORKS
#    def kandor_cooler(self):
#        if self.kandor_cooler_checkBox.isChecked():
#            self.myKandor.cooler = True
#            print('Cooler ON')
#        else:
#            self.myKandor.cooler = False
#            print('Cooler OFF')
#jks68 20/10/2021 end          
    def set_kandor_grating(self):
        self.myKandor.kymera.SetGrating(grating_num=int(self.TriaxGratingNumber_comboBox.currentText()))
        
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
        
    def laser633_isON(self):
        if self.laser633_isON_checkbox.isChecked():
            self.laser_633.turn_on()
        else:
            self.laser_633.turn_off()

    def update_Raman_spectrum_plot(self, spectrum):
        self.RamanWavelengths = self.myKandor.get_x_axis()
#        self.RamanWavelengths = self.KandorControlUI.Andor.x_axis       
        self.RamanSpectrum_plot.plot(self.RamanWavelengths, spectrum, clear = True, pen = 'r')  # plot current Raman spectrum in real time, ##sunny added 'np.flip(spectrum),' as the graph was flipped over the centralw wavelength
        self.RamanSpectrumImagePlotData.append(spectrum)     # plot Raman spectra over time as image plot        
        self.RamanImagePlotItem.setImage(np.asarray(self.RamanSpectrumImagePlotData))
               
    def update_darkfield_spectrum_plot(self, spectrum):
        self.darkfieldSpectrum_plot.plot(self.darkfieldWavelengths, spectrum, clear = True, pen = 'r')  # plot current darkfield spectrum in real time
        self.darkfieldSpectrumImagePlotData.append(np.nan_to_num(spectrum)[80:750])     # plot darkfield-time in range 400-900nm as image plot. ImageItem doesn't handle nan values, converting those to zero   
        self.darkfieldImagePlotItem.setImage(np.asarray(self.darkfieldSpectrumImagePlotData))
        
    def update_electronic_plot(self, timePlotInput, voltagePlotInput, currentPlotInput):
        if not hasattr(self,'voltages_data'):
            self.electronics_plot.clear()
            self.electronicsIV_plot.clear()
            self.times_data = [timePlotInput]
            self.voltages_data = [voltagePlotInput]
            self.currents_data = [currentPlotInput]
#            self.capacitance_data = [0]
        else:
            self.times_data.append(timePlotInput)
            self.voltages_data.append(voltagePlotInput)
            self.currents_data.append(currentPlotInput)
#            if (len(self.voltages_data)>1):### we added this line and the next 2
#                capacitance=currentPlotInput/((self.voltages_data[-1]-self.voltages_data[-2])/(self.times_data[-1]-self.times_data[-2]))
#                self.capacitance_data.append(capacitance)
#              #  print(str(capacitance))
#            else:
#                self.capacitance_data.append(0)
#                
#        self.capacitance_plot.plot(self.voltages_data, self.capacitance_data, clear = True, pen = 'r')
        self.electronics_plot.plot(self.times_data, self.currents_data, clear = True, pen = 'r')
        self.electronics_IVplot.plot(self.voltages_data, self.currents_data, clear = True, pen = 'r')

    def rampChangeDirection(self):
        self.voltageRampSign *= -1  
    
    def open_SmarAct_UI(self):
        delattr(self,'SmarAct_stage')
        self.SmarAct_stage = SmaractMCSSerial('COM6',3)
        self.SmarAct_stage.show_gui(blocking=False)

    def open_Andor_UI(self):
#        self.KandorControlUI = self.myKandor.show_gui(block = False)
        self.KandorControlUI = self.myKandor.get_control_widget()
        self.KandorPreviewUI = self.myKandor.get_preview_widget()
        self.KandorControlUI.show()
        self.KandorPreviewUI.show()

    def open_laser633_UI(self):
        self.laser633.show_gui()
        
    def open_SMC100_ui(self):
        self.SMC100.show_gui()

    def open_Dawn_z_stack_ui(self):
        self.z_stack.show_gui()
        
    def acquireIVdatapoint(self, activeVoltage, t0, activeDatagroup):
        measuredCurrent = self.smu.read_current()
        electronicActiveTime = (time.time()-t0)
        activeDatagroup.append_dataset("voltages", activeVoltage)
        activeDatagroup.append_dataset("currents", measuredCurrent)
        activeDatagroup.append_dataset("electronic_times", electronicActiveTime)
#        if len(self.capacitance_data)==0:
#            activeDatagroup.append_dataset("Capacitance", 0)
#        else:
#            activeDatagroup.append_dataset("Capacitance", self.capacitance_data[-1])
        self.live_electronic_signal.emit(electronicActiveTime, activeVoltage, measuredCurrent)
        
    
    def acquire_Raman_spectrum(self):
        self.RamanSpectrum = np.asarray(self.myKandor.capture()[0] )  # capture() retruns a tuple, converting to numpy array
        
    def acquire_darkfield_spectrum(self):
        self.darkfieldSpectrum = np.asarray(self.OOspectrometer.read_processed_spectrum())
        
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
    
    def flipShutter1(self):
        return

    def flipShutter2(self):
        return

    def flipMirror(self):
        return
    
    def acquire_single_Raman_spectrum(self):
        self.myKandor.set_camera_parameter('Exposure', self.RamanIntegrationTime)
        self.myKandor.AcquisitionMode = 1
        self.myKandor.ReadMode = 3
        self.myKandor.set_camera_parameter('SingleTrack', self.centre_row, self.num_rows)
        self.myKandor.kymera.SetWavelength(self.shamrockWavelength_spinBox.value())
        self.RamanWavelengths = self.myKandor.get_x_axis()
#        self.RamanWavelengths = self.KandorControlUI.Andor.x_axis
        self.RamanBackground = self.KandorControlUI.Andor.background
    
        time.sleep(0.5)
        self.RamanSpectrum = np.asarray( self.myKandor.capture()[0] )
        self.RamanSpectrum_plot.plot(self.RamanWavelengths, self.RamanSpectrum, clear = True, pen = 'r')

    def save_single_Raman_spectrum(self):
        activeSingleRamanDataset = self.singleRamanSpectraGroup.create_dataset('singleRamanSpectrum_%d', data = self.RamanSpectrum)
        activeSingleRamanDataset.attrs.create("Description", str(self.single_Raman_spectrum_description))
        activeSingleRamanDataset.attrs.create("Raman_Wavelengths", self.RamanWavelengths)
        activeSingleRamanDataset.attrs.create('Raman_Integration_Time', self.RamanIntegrationTime)
        activeSingleRamanDataset.attrs.create('Raman_Background', self.RamanBackground)

                
    def shutdown(self):
        self.activeDatafile.close()
        self.myKandor.cooler = False
        self.myKandor.close()
        self.OOspectrometer.shutdown_seabreeze()
        print('----Experiment ended----')
        #TODO - closes any UI windows
    
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
        


if __name__ == '__main__':
    activeDatafile = nplab.current_datafile(working_directory = os.path.abspath(os.path.join(os.getcwd(),"../../..")))
    #working directory should now be C:\\Users\\<name>\\Documents    
    gui_activeDatafile = activeDatafile.get_qt_ui()
    gui_activeDatafile.show()
    
    experiment = Lab3_experiment(activeDatafile)
    experiment.show_gui()
    print('Done')
