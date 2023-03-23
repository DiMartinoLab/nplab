"""
Dawn - dmk50
October 2022 - aiming to have new gui with control of both ocean optics and the stage - to enable df-stack
Trung - xtn20
November 2022 - SmarAct stage control with MCS2 controller
"""
#nb, next up - look at spectrometer displayui

from PyQt5 import QtWidgets, uic
# enable the parallel processing
from PyQt5.QtCore import QRunnable, Qt, QThreadPool
import sys
import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from nplab.utils.gui import QtCore, QtGui, uic, get_qt_app, show_widget
from nplab.instrument.spectrometer import Spectrometer, SpectrometerDisplayUI
from nplab.experiment import Experiment, ExperimentStopped
from nplab.instrument.spectrometer.seabreeze import OceanOpticsSpectrometer, OceanOpticsControlUI
import smaract.ctl as smaract_package
from nplab.instrument.stage.SMC100 import SMC100 as SMC_package
from nplab.ui.ui_tools import UiTools
import ctypes
from ctypes import byref, c_int, c_ulong, c_double
import pyqtgraph as pg
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg
from matplotlib.figure import Figure
from PySide6.QtCharts import QChart, QChartView, QLineSeries, QDateTimeAxis, QValueAxis
import time

#import curve to be used when using reflective metal as reference material rather than white scatterer
ratio_curve_points = np.genfromtxt(fname='ratio_curve.txt')

class z_stack_window_object(Experiment, QtWidgets.QWidget, UiTools):
    #inherit QWidget methods
    def __init__(self, activeDatafile, OOSpect_instance, Stage_instance, MiTC_instance, MiPS_instance, stage_enabled,  ui_file = os.path.join(os.path.dirname(__file__),'z_stack_window.ui'),   parent=None):
        
        # set objects for spectrometer and stage
        self.OOspectrometer = OOSpect_instance
        self.OOspectrometer.set_integration_time(50)
        self.stage_enabled = stage_enabled
        if self.stage_enabled == 'Cryostat stage':
            self.Smaract_stage = Stage_instance
        elif self.stage_enabled == 'Olympus stage':
            self.SMC_stage = Stage_instance

        # set objects for the MiTC and MiPS
        self.MiTC_handle = MiTC_instance
        self.MiPS_handle = MiPS_instance
        
        # variable for checking if communication is free or busy
        global MiTC_communication
        MiTC_communication = 0
        global MiPS_communication
        MiPS_communication = 0
        # devices
        if self.MiTC_handle != None and self.MiPS_handle != None : 
            self.MiTC_handle.devices = {'Magnet' : 'DEV:DB1.T1:TEMP', 'HTS' : 'DEV:DB2.T1:TEMP', 'Sample' : 'DEV:MB1.T1:TEMP', 'Heater' : 'DEV:MB0.H1:HTR'}
            self.MiPS_handle.devices = {'Magnet' : 'DEV:GRPZ:PSU'}
            self.MiPS_handle.getMiPSHeater = 'OFF'
        
        # make this window available for the main program
        super(z_stack_window_object, self).__init__() 
        uic.loadUi(ui_file, self)
        
        #name a group within the h5 file to save collected data into
        self.StackGUI_df_group = activeDatafile.create_group(name='Z-Stack_DataGroup', auto_increment=True)
        self.StackGUI_df_group.attrs.create('Version Date', str('Date of latest edit of Stack Gui = 21-11-2022'))
        #save spectra here, with addintional info such as stage location and cryostat params also saved.
        self.OO_Spectra_group= self.StackGUI_df_group.create_group(name='OO Spectra', auto_increment=True)
        self.Stacked_Spectra_group = self.StackGUI_df_group.create_group(name='Stacked Spectra', auto_increment=True)
        self.ratio_curve_values = [item[1] for item in ratio_curve_points]
        self.OO_Spectra_group.attrs.create('Ratio_curve', self.ratio_curve_values)
        #save cryostat temperature vs time data in self.Cryo_Temps
        self.Cryo_Temps= self.StackGUI_df_group.create_group('Cryostat Temperatures')
        """TO DO NOTE=  for all OOspectrometer quantities - match object attrubute naming conventions, 
        so that bkg and ref spectra will be used for both OO gui and z-stack gui"""
        
    
        """connect control and readout buttons for spectrometer temperature & spectrometer integration time """
        self.set_tec_temperature_pushButton.clicked.connect(self.gui_set_tec_temperature)
        self.read_tec_temperature_pushButton.clicked.connect(self.gui_read_tec_tempeature)
        initial_temperature = np.round(self.OOspectrometer.get_tec_temperature(), decimals = 1)
        self.tec_temperature_lcdNumber.display(float(initial_temperature))
        self.set_tec_temperature_LineEdit.setText(str(initial_temperature))
        self.set_integration_time_push_button.clicked.connect(self.gui_set_integration_time)
        
        """setup spectrometer commands """
        self.wl_in_nm = self.OOspectrometer.get_wavelengths()
        self.referencing_type = 'White'
        self.read_bkg_for_ref_button.clicked.connect(self.read_bkg_for_ref_method)
        self.clear_bkg_for_ref_button.clicked.connect(self.clear_bkg_for_ref_method)
        self.take_bkg_button.clicked.connect(self.read_background_method)
        self.clear_bkg_button.clicked.connect(self.clear_background_method)
        self.read_ref_button.clicked.connect(self.read_ref_method)
        self.clear_ref_button.clicked.connect(self.clear_ref_method)
        self.gold_ref_checkbox.stateChanged.connect(self.checkbox_testing)
        self.take_spectrum_button.clicked.connect(self.take_spec_and_update_display)

        
        self.save_button.clicked.connect(self.save_single_DF_spectrum)
        #self.averaging_enabled=False
        
        
        
        """setup display for spectrum plot """
        self.plotbox = QtWidgets.QGroupBox()
        self.plotbox.setLayout(QtWidgets.QGridLayout())
        self.plotlayout = self.plotbox.layout()     
        self.spec_plot = pg.PlotWidget(labels = {'bottom':'Wavelength (nm)'})
        self.plotlayout.addWidget(self.spec_plot)
        self.figure_widget = self.replace_widget(self.display_layout, self.figure_widget, self.plotbox) 
        
        """setup z-stack controls and display """
        self.take_z_stack_button.clicked.connect(self.z_stack_method)
        self.plotbox_stack  = QtWidgets.QGroupBox()
        self.plotbox_stack.setLayout(QtWidgets.QGridLayout())
        self.plotlayout_stack  = self.plotbox_stack.layout()     
        self.stack_plot=pg.PlotWidget(labels = {'bottom':'Wavelength (nm)'})
        self.plotlayout_stack .addWidget(self.stack_plot)
        self.figure_widget_stack = self.replace_widget(self.display_layout, self.figure_widget_stack, self.plotbox_stack ) 
        
        
        # Trung Nov. 2022
        
        # Check which stage is loaded 'Smaract at Cryostat' or 'SMC at Olympus'
        if stage_enabled == 'Cryostat stage':
            # added parts for the SmarAct stage control
            # Further infos and code: C:\SmarAct\MCS2\SDK\Python\examples
            # get information about the stage, number of channels
            self.smaract_serial = smaract_package.GetProperty_s(self.Smaract_stage, 0, smaract_package.Property.DEVICE_SERIAL_NUMBER)
            print("Device Serial Number: {}".format(self.smaract_serial))
            self.smaract_name = smaract_package.GetProperty_s(self.Smaract_stage, 0, smaract_package.Property.DEVICE_NAME)
            print("Device Name: {}".format(self.smaract_name))
            
            # Read channel info
            self.smaract_no_of_bus_modules = smaract_package.GetProperty_i32(self.Smaract_stage, 0, smaract_package.Property.NUMBER_OF_BUS_MODULES)
            print("Number of Bus Modules: {}".format(self.smaract_no_of_bus_modules))
            self.smaract_no_of_channels = smaract_package.GetProperty_i32(self.Smaract_stage, 0, smaract_package.Property.NUMBER_OF_CHANNELS)
            print("Total Number of Channels: {}".format(self.smaract_no_of_channels))
            # The "channelIndex" must be passed to all module properties. Note that it is zero-based.
            for i in range(self.smaract_no_of_channels):
                print("        Channel: {}".format(i))
                self.smaract_name = smaract_package.GetProperty_s(self.Smaract_stage, i, smaract_package.Property.POSITIONER_TYPE_NAME)  
                self.smaract_pos_type = smaract_package.GetProperty_i32(self.Smaract_stage, i, smaract_package.Property.POSITIONER_TYPE)
                print("        Positioner Type: {} ({})".format(self.smaract_name, self.smaract_pos_type))
                self.smaract_state = smaract_package.GetProperty_i32(self.Smaract_stage, i, smaract_package.Property.CHANNEL_STATE)  
                print("        Amplifier enabled: ", end='')
                if (self.smaract_state & smaract_package.ChannelState.AMPLIFIER_ENABLED) != 0:
                    print("yes")
                else:
                    print("no")
                # Here we read the "Channel Type" property to determine which additional properties and flags are of interest.
                # E.g. the maxCLF property is only available for Stick-Slip-Driver channels and the "isPhased"
                # Channel State flag is only available for Magnetic-Driver channels.
                # Note that the "Module Type" and "Channel Type" properties share the same list of types.
                self.smaract_type = smaract_package.GetProperty_i32(self.Smaract_stage, i, smaract_package.Property.CHANNEL_TYPE)
                if self.smaract_type == smaract_package.ChannelModuleType.STICK_SLIP_PIEZO_DRIVER:
                    self.smaract_max_clf = smaract_package.GetProperty_i32(self.Smaract_stage, i, smaract_package.Property.MAX_CL_FREQUENCY)
                    print("        Max-CLF: {} Hz".format(self.smaract_max_clf))
                elif self.smaract_type == smaract_package.ChannelModuleType.MAGNETIC_DRIVER:
                    print("        Channel is phased: ", end='')
                    if (self.smaract_state & smaract_package.ChannelState.IS_PHASED) != 0:
                        print("yes")
                    else:
                        print("no")
                print("-------------------------------------------------------")
            # set movement mode of the stage to open-loop STEP
            self.smaract_movemode = smaract_package.MoveMode.STEP
            # Set move mode depending properties for the next movement.
            # move_mode == ctl.MoveMode.STEP:
            # Set step frequency [in Hz].
            # Valid range: 1 to 20000 Hz
            
            # Set maximum step amplitude [in dac increments].
            # valid range: 0 to 65535 corresponding to 0 to 100V piezo voltage
            # Lower amplitude values result in smaller step width.
            
            # set for all 3 channels
            # channel 0/2 for left-right/up-down depending on which point of view
            # see "current_view" for more infos
            channel = 0
            self.smaract_set_movemode(channel, self.smaract_movemode)
            channel = 1
            self.smaract_set_movemode(channel, self.smaract_movemode)
            channel = 2
            self.smaract_set_movemode(channel, self.smaract_movemode)
            
            # set view
            self.stage_current_view_button.setText('Camera view')
            
            self.smaract_channel_up = 0
            self.smaract_direction_up = -1
            
            self.smaract_channel_down = 0
            self.smaract_direction_down = 1
            
            self.smaract_channel_left = 2
            self.smaract_direction_left = -1
            
            self.smaract_channel_right = 2
            self.smaract_direction_right = 1
    
            
            # END stage property
            
        
            # link buttons for smaract movement to functions
            self.stage_mleft_button.clicked.connect(self.smaract_moveleft)
            self.stage_mright_button.clicked.connect(self.smaract_moveright)
            self.stage_mup_button.clicked.connect(self.smaract_moveup)
            self.stage_mdown_button.clicked.connect(self.smaract_movedown)
            
            self.stage_mtoward_button.clicked.connect(self.smaract_movetoward)
            self.stage_maway_button.clicked.connect(self.smaract_moveaway)
        
            # Movement of stage on cryostat may differ from movement on camera view
            # Button to change view corresponding to the movement name
            self.stage_current_view_button.clicked.connect(self.current_view)
        
        # Disconnect Smaract and Ocean Optics
        self.stageOO_disconnect_button.clicked.connect(self.stageOO_disconnect)
        
        # Mercury iTC and iPS-M
        self.SetSampleTemp_button.clicked.connect(self.SetSampleTemp)
        self.MercuryiPSHeater_button.clicked.connect(self.MercuryiPSHeater)
        self.setMagnetSetPoint_button.clicked.connect(self.setMagnetSetPoint)
        self.setMagnetRamp_button.clicked.connect(self.setMagnetRamp)
        self.MagnetRead_button.clicked.connect(self.MagnetRead)
        self.MagnetGo_button.clicked.connect(self.MagnetGo)
        self.MagnetToZero_button.clicked.connect(self.MagnetToZero)
        self.ReadTemperatures_button.clicked.connect(self.ReadTemperatures)
        self.Hold_button.clicked.connect(self.MagnetHold)
        
        """ parallel processing"""
        self.threadCount = QThreadPool.globalInstance().maxThreadCount()
        self.pool = QThreadPool.globalInstance()
        
        
    def read_bkg_for_ref_method(self):
        #read spec and store as background for reference. chenge check-box state to match.
        #self.OOspectrometer.read_background_ref
        self.OOspectrometer.background_ref = self.OOspectrometer.read_spectrum()
        self.OOspectrometer.background_for_ref_int_t = self.OOspectrometer.integration_time
        self.OOspectrometer.stored_backgrounds[self.OOspectrometer.reference_ID] = {'background_ref' : self.OOspectrometer.background_ref,'background_for_ref_int_t': self.OOspectrometer.background_for_ref_int_t}
        self.OOspectrometer.update_config('background_ref', self.OOspectrometer.background_ref)
        self.OOspectrometer.update_config('background_for_ref_int_t', self.OOspectrometer.background_for_ref_int_t)
        print(self.OOspectrometer.read_background_ref)
        self.background_subtracted_ref.setCheckState(QtCore.Qt.Checked)
            
    def clear_bkg_for_ref_method(self):
        self.OOspectrometer.clear_background_ref
        self.background_subtracted_ref.setCheckState(QtCore.Qt.Unchecked)
            
    def read_ref_method(self):
        #self.OOspectrometer.read_reference
        self.local_reference = self.OOspectrometer.read_spectrum() 
        self.OOspectrometer.reference = self.local_reference
        self.OOspectrometer.reference_int_t = self.OOspectrometer.integration_time
        self.OOspectrometer.update_config('reference', self.OOspectrometer.reference)
        self.OOspectrometer.update_config('reference_int',self.OOspectrometer.reference_int_t ) 
        self.OOspectrometer.stored_references[self.OOspectrometer.reference_ID] = {'reference' : self.OOspectrometer.reference,'reference_int_t' : self.OOspectrometer.reference_int_t}
        self.referenced.setCheckState(QtCore.Qt.Checked)

    def clear_ref_method(self):
        self.OOspectrometer.clear_reference
        self.referenced.setCheckState(QtCore.Qt.Unchecked)

    def read_background_method(self):
        """Acquire a new spectrum and use it as a background measurement.
        This background should be less than 50% of the spectrometer saturation"""
        self.OOspectrometer.background = self.OOspectrometer.read_spectrum()
        self.OOspectrometer.background_int_t = self.OOspectrometer.integration_time
        self.OOspectrometer.stored_backgrounds[self.OOspectrometer.reference_ID] = {'background' : self.OOspectrometer.background,
                                                     'background_int': self.OOspectrometer.background_int_t}
        self.OOspectrometer.update_config('background', self.OOspectrometer.background)
        self.OOspectrometer.update_config('background_int', self.OOspectrometer.background_int_t)
        self.background_subtracted.setCheckState(QtCore.Qt.Checked)
        
        
    def clear_background_method(self):
        self.OOspectrometer.background = None
        self.OOspectrometer.background_int_t = None
        self.background_subtracted.setCheckState(QtCore.Qt.Unchecked)
    
    def take_spec(self):

        #update int t
        print('take spec using int t = '+str(self.OOspectrometer.integration_time))
        # parallel processing - variable for checking if taking spectrum is done
        global finished_OO
        finished_OO = 0
        # read the current field strength and sample temperature
        if self.MiPS_handle != None:
            current_field = sub_read_current_field(self.MiPS_handle, self.cFieldT)
            self.pool.start(current_field)
        if self.MiTC_handle  != None:
            current_temp = sub_read_current_sample_temp(self.MiTC_handle, self.cSampleTemp)
            self.pool.start(current_temp)
        
        if self.stage_enabled == 'Cryostat stage':
            # change the movemode of channel 1 to absolute and
            # read the current stage position (focussing channel: 1)
            channel = 1
            move_mode = smaract_package.MoveMode.CL_ABSOLUTE
            self.smaract_set_movemode(1, move_mode)
            self.current_z_position = smaract_package.GetProperty_i64(self.Smaract_stage, channel, smaract_package.Property.POSITION)
            
#            # test
#            smaract_package.Move(self.Smaract_stage, channel, 0, 0)
#            time.sleep(5)
#            #end test
            
            # hold the stage at this position while OO takes spectrum
            stage_on_hold = sub_stage_on_hold(self.Smaract_stage, self.current_z_position , channel)
            self.pool.start(stage_on_hold)
            
        # Take the spectrum in the background
        take_spec_task = sub_take_spec(self.OOspectrometer)
        self.pool.start(take_spec_task)
        
        # Wait till spectrum is taken
        while finished_OO < 1:
            just_a_random_variable = 0
            
        if self.stage_enabled == 'Cryostat stage':
            # change the movemode of channel 1 back to step
            channel = 1
            self.smaract_set_movemode(channel, self.smaract_movemode)
        
        self.current_raw_spec = take_spec_task.current_spec
        if self.MiPS_handle != None:
            self.current_field = current_field.cFieldT_value
        elif self.MiPS_handle == None:
            self.current_field =0
        if self.MiTC_handle != None:  
            self.current_temp = current_temp.cSampleTemp_value
        elif self.MiTC_handle == None: 
            self.current_temp  = 0
            
        print(str(self.current_field) + ' T and ' + str(self.current_temp) + ' K')        
        #storing raw and processed spec in seperate variables
        self.current_processed_spec = self.process_spec(self.current_raw_spec)
        
    def take_spec_and_update_display(self):
        # I suggest to use "take_spec" function here instead of "read_spectrum"
#        self.current_spec = self.OOspectrometer.read_spectrum()
        self.take_spec()
        #this line in twice to get rid of weird timing glitch...temporary patch...
        self.take_spec()
        global finished_OO
        while finished_OO < 1:
            test = 0
        if finished_OO==1:
            self.update_display(self.current_processed_spec)
        
    def process_spec(self, a_spectrum):
        """check if bkg for ref exists. if so, scale to match spec int time. if not, use placeholder array of zeros."""
        if np.array(self.OOspectrometer.background_ref).all() != None:
            scaled_bkg_for_ref = self.OOspectrometer.background_ref * (self.OOspectrometer.integration_time/self.OOspectrometer.background_for_ref_int_t)
        elif self.OOspectrometer.background_ref ==None:
            scaled_bkg_for_ref  = np.zeros(len(a_spectrum))
        """check if bkg  exists. if so, scale to match spec int time. if not, use placeholder array of zeros."""
        if np.array(self.OOspectrometer.background).all() != None:
            scaled_bkg = self.OOspectrometer.background * (self.OOspectrometer.integration_time/self.OOspectrometer.background_int_t)
        elif self.OOspectrometer.background == None:
            scaled_bkg = np.zeros(len(a_spectrum))
        """check if reference exists. if so, scale to match spec int time and multiply by ratio curve if "Au ref" box is ticked. """
        """if not, replace with array of ones - then calculated processed spectrum"""
        if self.referencing_type == 'Metal' and np.array(self.OOspectrometer.reference).all() != None:
            scaled_ref =self.OOspectrometer.reference * self.ratio_curve_values * (self.OOspectrometer.integration_time/self.OOspectrometer.reference_int_t)
            processed_spectrum = (a_spectrum - scaled_bkg)/(scaled_ref - scaled_bkg_for_ref)
        elif self.referencing_type == 'White' and np.array(self.OOspectrometer.reference).all() != None:
            scaled_ref =self.OOspectrometer.reference * (self.OOspectrometer.integration_time/self.OOspectrometer.reference_int_t)
            processed_spectrum = (a_spectrum - scaled_bkg)/(self.OOspectrometer.reference- scaled_bkg_for_ref)
        elif self.OOspectrometer.reference == None:
            ref_placeholder = np.ones(len(a_spectrum))
            processed_spectrum = (a_spectrum - scaled_bkg)/(ref_placeholder- scaled_bkg_for_ref)
        return processed_spectrum

        
    def button_test(self):
            print('button is connected')
          
    def update_display(self, spectrum):
        #replace any Nan in input spectrum with a 0
        spectrum = np.array([0 if np.isnan(i) else i for i in spectrum])
        #empty the display each time so that only most recent plot is shown
        self.spec_plot.clear()
        self.spec_plot.plot(self.wl_in_nm, spectrum)

    """connect spectrometer temperature control and readout buttons"""
    def gui_set_tec_temperature(self):
        self.OOspectrometer.tec_temperature = float(self.set_tec_temperature_LineEdit.text().strip())
    def gui_read_tec_tempeature(self):
        self.tec_temperature_lcdNumber.display(float(self.OOspectrometer.tec_temperature)) 
    """connect spectrometer integration time control and readout buttons""" 
    def gui_set_integration_time(self):
        self.OOspectrometer.integration_time = float(self.set_integration_time_LineEdit.text().strip())
        self.OOspectrometer.set_integration_time(self.OOspectrometer.integration_time)
        print('Intergration time set to: '+str(self.OOspectrometer.integration_time))
        
        
    
    
    def checkbox_testing(self, state):
        sender = self.sender()
        if sender is self.gold_ref_checkbox and state == QtCore.Qt.Checked:
            print('gold ref checkbox is checked')
            self.referencing_type = 'Metal'
        elif sender is self.gold_ref_checkbox and state == QtCore.Qt.Unchecked:
            print('gold ref checkbox is NOT checked')
            self.referencing_type = 'White'
        elif sender is self.average_checkBox and state == QtCore.Qt.Checked:
            self.OOspectrometer.averaging_enabled =True
        elif sender is self.average_checkBox and state == QtCore.Qt.Unchecked:
            self.OOspectrometer.averaging_enabled =False
            
    def read_spec_description(self):
            self.current_spec_description= str(self.description_LineEdit.text().strip())
            
    def z_stack_method(self):
        self.stack_plot.clear()
        self.stack_step_size = int(self.zstep_box_2.value())
        self.stack_num_steps = int(self.num_steps_box.value())
        initial_z = smaract_package.GetProperty_i64(self.Smaract_stage, 1, smaract_package.Property.POSITION)
        global z_stack_running
        z_stack_running = 1
        step_num=0
        self.processed_stack_spectra = []
        self.raw_stack_spectra = []
        while step_num<self.stack_num_steps:
            self.take_spec()
            global finished_OO
            while finished_OO < 1:
                test = 0
            if finished_OO==1:
            #self.update_display(self.current_processed_spec)
            #self.current_raw_spec = take_spec_task.current_spec
            # Holding stage: ensure the stage is not moving during measurement?? TO DO
                self.smaract_move(1, self.stack_step_size*(-1))
                self.processed_stack_spectra = self.processed_stack_spectra + [self.current_processed_spec]
                self.raw_stack_spectra = self.raw_stack_spectra + [self.current_raw_spec]
                self.stack_plot.plot(self.wl_in_nm, self.current_processed_spec)
                step_num = step_num+1
            
        activeStackDataset = self.Stacked_Spectra_group.create_dataset('stackedDFSpectra_%d', data = self.raw_stack_spectra)
        activeStackDataset.attrs.create('stack step size', self.stack_step_size)
        activeStackDataset.attrs.create('stack step number', self.stack_num_steps)
        activeStackDataset.attrs.create("Wavelengths", self.wl_in_nm)
        activeStackDataset.attrs.create('Integration_Time', self.OOspectrometer.integration_time)
        activeStackDataset.attrs.create("Initial z position", initial_z)
        final_z = smaract_package.GetProperty_i64(self.Smaract_stage, 1, smaract_package.Property.POSITION)
        activeStackDataset.attrs.create("Final z position", final_z)
        self.read_spec_description()
        activeStackDataset.attrs.create("Description", self.current_spec_description)
        activeStackDataset.attrs.create("Processed Spectra", self.processed_stack_spectra )
        
        activeStackDataset.attrs.create('bkg_for_ref_int_t', self.OOspectrometer.background_for_ref_int_t)
        activeStackDataset.attrs.create('bkg_int_t', self.OOspectrometer.background_int_t)
        activeStackDataset.attrs.create('ref_int_t', self.OOspectrometer.reference_int_t)
        activeStackDataset.attrs.create('bkg_for_ref', self.OOspectrometer.background_ref )
        activeStackDataset.attrs.create('background', self.OOspectrometer.background)
        activeStackDataset.attrs.create('reference', self.OOspectrometer.reference)
        
        activeStackDataset.attrs.create('solenoid temperature', self.current_temp)
        activeStackDataset.attrs.create('solenoid field', self.current_field)
        

        #TRUNG TO DO TEST IPS
        
        
        self.make_stack_plot()
            
    def make_stack_plot(self):
        spec_nums = np.arange(0,self.stack_num_steps,1)
        spectra_T = np.transpose(self.processed_stack_spectra)
        plt.pcolormesh(spec_nums, self.wl_in_nm, spectra_T, shading = 'auto', norm = colors.SymLogNorm(base = np.e, linthresh = 0.5, linscale = 0.05, vmin = -0.5, vmax = 0.5))
        plt.clim(0,0.1)
        plt.show()
        
    
    def save_single_DF_spectrum(self):
        #save a raw spectrum = self.current_raw_spec
        """nb - want to save both raw and processed!"""
        activeSingleDFDataset = self.OO_Spectra_group.create_dataset('singleDFSpectrum_%d', data = self.current_raw_spec)

        self.read_spec_description()
        activeSingleDFDataset.attrs.create("Description", self.current_spec_description)
        activeSingleDFDataset.attrs.create("Processed Spectrum", self.current_processed_spec)
        activeSingleDFDataset.attrs.create("Wavelengths", self.wl_in_nm)
        activeSingleDFDataset.attrs.create('Integration_Time', self.OOspectrometer.integration_time)
    
        activeSingleDFDataset.attrs.create('bkg_for_ref_int_t', self.OOspectrometer.background_for_ref_int_t)
        activeSingleDFDataset.attrs.create('bkg_int_t', self.OOspectrometer.background_int_t)
        activeSingleDFDataset.attrs.create('ref_int_t', self.OOspectrometer.reference_int_t)
        
        activeSingleDFDataset.attrs.create('bkg_for_ref', self.OOspectrometer.background_ref )
        activeSingleDFDataset.attrs.create('background', self.OOspectrometer.background)
        activeSingleDFDataset.attrs.create('reference', self.OOspectrometer.reference)
        
        activeSingleDFDataset.attrs.create('solenoid temperature', self.current_temp)
        activeSingleDFDataset.attrs.create('solenoid field', self.current_field)
        
        activeSingleDFDataset.attrs.create('stage z position', self.current_z_position )
        
        
        
    # Trung Nov. 2022
    # added parts for the SmarAct stage control
    # Further infos and code: C:\SmarAct\MCS2\SDK\Python\examples
    # self.Smaract_stage = d_handle
  
    
    def smaract_set_movemode(self, channel, smaract_movemode):
        if smaract_movemode == smaract_package.MoveMode.STEP:
            # Change movement mode to step movement
            # Set maximum step amplitude [in dac increments].
            # valid range: 0 to 65535 corresponding to 0 to 100V piezo voltage
            # Lower amplitude values result in smaller step width.
            smaract_package.SetProperty_i32(self.Smaract_stage, channel, smaract_package.Property.MOVE_MODE, smaract_movemode)
            smaract_package.SetProperty_i32(self.Smaract_stage, channel, smaract_package.Property.STEP_FREQUENCY, 1000)
            smaract_package.SetProperty_i32(self.Smaract_stage, channel, smaract_package.Property.STEP_AMPLITUDE, 65535)
        if smaract_movemode == smaract_package.MoveMode.CL_ABSOLUTE:
            # Change movement mode to absolute
            smaract_package.SetProperty_i32(self.Smaract_stage, channel, smaract_package.Property.MOVE_MODE, smaract_movemode)
            # Set move velocity [in pm/s].
            smaract_package.SetProperty_i64(self.Smaract_stage, channel, smaract_package.Property.MOVE_VELOCITY, 1000000000)
            # Set move acceleration [in pm/s2].
            smaract_package.SetProperty_i64(self.Smaract_stage, channel, smaract_package.Property.MOVE_ACCELERATION, 1000000000)
                  
    # move a channel
    def smaract_move(self, channel, move_value):
        # Specify the number of steps to perform and the direction.
        print("MCS2 open loop step move, channel {}, steps: {}.".format(channel, move_value))
        # Start actual movement.
        smaract_package.Move(self.Smaract_stage, channel, move_value, 0)
        
    # Note that the function call returns immediately, without waiting for the movement to complete.
    # The "ChannelState.ACTIVELY_MOVING" (and "ChannelState.CLOSED_LOOP_ACTIVE") flag in the channel state
    # can be monitored to determine the end of the movement.
    # stop a channel
    def stop(self, channel):
        print("MCS2 stop channel: {}.".format(channel))
        smaract_package.Stop(self.Smaract_stage, channel)

    def smaract_moveleft(self):
        # channel 0 (+)
        self.smaract_move_value = int(self.xystep_Box.value())
        self.smaract_move(self.smaract_channel_left, self.smaract_move_value*self.smaract_direction_left)
        
    def smaract_moveright(self):
        # channel 0 (-)
        self.smaract_move_value = int(self.xystep_Box.value())
        self.smaract_move(self.smaract_channel_right, self.smaract_move_value*self.smaract_direction_right)
        
    def smaract_moveup(self):
        # channel 2 (-)
        self.smaract_move_value = int(self.xystep_Box.value())
        self.smaract_move(self.smaract_channel_up, self.smaract_move_value*self.smaract_direction_up)
        
    def smaract_movedown(self):
        # channel 2 (+)
        self.smaract_move_value = int(self.xystep_Box.value())
        self.smaract_move(self.smaract_channel_down, self.smaract_move_value*self.smaract_direction_down)
        
    def smaract_movetoward(self):
        # channel 1 (+)
        self.smaract_move_value = int(self.zstep_box.value())
        self.smaract_move(1, self.smaract_move_value)
        
    def smaract_moveaway(self):
        # channel 1 (-)
        self.smaract_move_value = int(self.zstep_box.value())
        self.smaract_move(1, self.smaract_move_value*(-1))
    
    # Switch stage channel corresponding to the view
    def current_view(self):
        button_text = self.stage_current_view_button.text()
        if button_text == 'Camera view':
            self.stage_current_view_button.setText('Cryostat view')
            
            self.smaract_channel_left = 0
            self.smaract_direction_left = 1
            
            self.smaract_channel_right = 0
            self.smaract_direction_right = -1
            
            self.smaract_channel_up = 2
            self.smaract_direction_up = -1
            
            self.smaract_channel_down = 2
            self.smaract_direction_down = 1
        
        elif button_text == 'Cryostat view':
            self.stage_current_view_button.setText('Camera view')
            
            self.smaract_channel_up = 0
            self.smaract_direction_up = 1
            
            self.smaract_channel_down = 0
            self.smaract_direction_down = -1
            
            self.smaract_channel_left = 2
            self.smaract_direction_left = 1
            
            self.smaract_channel_right = 2
            self.smaract_direction_right = -1
            
    
    def smaract_position(self):
        self.z_position = smaract_package.GetProperty_i64(self.Smaract_stage, 1, smaract_package.Property.POSITION)
        print("MCS2 position: {} pm.".format(self.z_position))
    
    
    # added parts for the SMC stage control
    



        
        
    # Control part for Mercury iTC and iPS
    
    def ReadTemperatures(self):
        background_event = sub_read_current_temps(self.MiTC_handle, self.cSampleTemp, self.cMagnetTemp, self.cHTSTemp)
        self.pool.start(background_event)
        
    def SetSampleTemp(self):
        background_event = sub_set_sample_temp(self.MiTC_handle, self.SetSampleTemp_box)
        self.pool.start(background_event)
        
    def MercuryiPSHeater(self):
        background_event = sub_heater_MiPS(self.MiPS_handle, self.MercuryiPSHeater_button)
        self.pool.start(background_event)
        
    def setMagnetSetPoint(self):
        background_event = sub_set_field(self.MiPS_handle, self.setMagnetSetPoint_box)
        self.pool.start(background_event)
        
    def setMagnetRamp(self):
        background_event = sub_set_ramp(self.MiPS_handle, self.setMagnetRamp_box)
        self.pool.start(background_event)
        
    def MagnetRead(self):
        background_event = sub_read_current_field(self.MiPS_handle, self.cFieldT)
        self.pool.start(background_event)
        
    def MagnetGo(self):
        background_event = sub_MiPS_GoToSetPoint(self.MiPS_handle)
        self.pool.start(background_event)
        
    def MagnetToZero(self):
        background_event = sub_MiPS_ToZero(self.MiPS_handle)
        self.pool.start(background_event)
    
    def MagnetHold(self):
        background_event = sub_MagnetHold(self.MiPS_handle)
        self.pool.start(background_event)
        
    def stageOO_disconnect(self):
        if self.stage_enabled == 'Cryostat stage':
            smaract_package.Close(self.Smaract_stage)
        elif self.stage_enabled == 'Olympus stage':
            print('nothing yet')
        self.OOspectrometer._close()
        print('Stage and Ocean Optics successfully disconnected')
    # END: added parts for the SmarAct stage control

    def make_window(self):
        app = get_qt_app()
        self.show()
        app.exec_()
        return self
    
    
    
    
    
    
    
    
    
    
    


class sub_take_spec(QRunnable):
    """this class uses a seperate thread while taking a spectrum."""
    def __init__(self, OOspectrometer):
        super().__init__()
        # this is local access to OOspectrometer. "self" here is the "sub_take_spec"
        self.OOspectrometer = OOspectrometer
    
    def run(self):
        self.current_spec = self.OOspectrometer.read_spectrum()
        #print(self.OOspectrometer.integration_time)
        # check if spectra taken is done to terminate stage correction
        global finished_OO
        finished_OO = 1
        print('done')
    
class sub_stage_on_hold(QRunnable):
    """this class uses a seperate thread while checking/controlling stage position. Runs while sub_take_spec is running
    (holds stage steady while spectrum is being collected) - stops once spectrum is finished"""
    def __init__(self, Smaract_stage, position, channel):
        super().__init__()
        # this is local access to Smaract_stage. "self" here is the "sub_stage_on_hold"
        self.Smaract_stage = Smaract_stage
        self.channel = channel
        self.position = position
    
    def run(self):
        global finished_OO
        while finished_OO < 1:
            cposition = smaract_package.GetProperty_i64(self.Smaract_stage, self.channel, smaract_package.Property.POSITION)
            if np.abs(cposition - self.position) > 500000: # = 500nm
                smaract_package.Move(self.Smaract_stage, self.channel, self.position, 0)
                print('Stage position corrected')

class sub_read_current_temps(QRunnable):
    def __init__(self, MiTC_handle, cSampleTempDisplay, cMagnetTempDisplay, cHTSTempDisplay):
        super().__init__()
        self.MiTC_handle = MiTC_handle
        # 3 display fields on the UI for sample, magnet and HTS
        self.cSampleTempDisplay = cSampleTempDisplay
        self.cMagnetTempDisplay = cMagnetTempDisplay
        self.cHTSTempDisplay = cHTSTempDisplay
    
    def run(self):
        # check and block the communication with MiTC while reading
        global MiTC_communication
        # check the communication is free, if not, wait till it is freed
        while MiTC_communication > 0:
            time.sleep(0.05)
        # block the communication
        MiTC_communication = 1
        
        # read temperature here
        self.cSampleTemp_value = str(self.MiTC_handle.getTemp('Sample'))
        self.cSampleTempDisplay.display(float(self.cSampleTemp_value))
        self.cMagnetTemp_value = str(self.MiTC_handle.getTemp('Magnet'))
        self.cMagnetTempDisplay.display(float(self.cMagnetTemp_value))
        self.cHTSTemp_value = str(self.MiTC_handle.getTemp('HTS'))
        self.cHTSTempDisplay.display(float(self.cHTSTemp_value))
        # free the communication
        MiTC_communication = 0
        print('Read temp done!')
        
class sub_read_current_sample_temp(QRunnable):
    def __init__(self, MiTC_handle, cSampleTempDisplay):
        super().__init__()
        self.MiTC_handle = MiTC_handle
        # 3 display fields on the UI for sample, magnet and HTS
        self.cSampleTempDisplay = cSampleTempDisplay
        self.cSampleTemp_value = self.cSampleTempDisplay.value()
    
    def run(self):
        # check and block the communication with MiTC while reading
        global MiTC_communication
        # check the communication is free, if not, wait till it is freed
        while MiTC_communication > 0:
            time.sleep(0.05)
        # block the communication
        MiTC_communication = 1
        # read temperature here
        self.cSampleTemp_value = str(self.MiTC_handle.getTemp('Sample'))
        self.cSampleTempDisplay.display(float(self.cSampleTemp_value))
        # free the communication
        MiTC_communication = 0
    
class sub_set_sample_temp(QRunnable):
    def __init__(self, MiTC_handle, SetSampleTemp_box):
        super().__init__()
        self.MiTC_handle = MiTC_handle
        self.setSampleTemp = float(SetSampleTemp_box.value())
    
    def run(self):
        # check and block the communication with MiTC while reading
        global MiTC_communication
        # check the communication is free, if not, wait till it is freed
        while MiTC_communication > 0:
            time.sleep(0.05)
        # block the communication
        MiTC_communication = 1
        # set temperature here
        """not yet doneeeeeeee"""
        # free the communication
        MiTC_communication = 0
        
class sub_heater_MiPS(QRunnable):
    def __init__(self, MiPS_handle, MiPSHeater_button):
        super().__init__()
        self.MiPS_handle = MiPS_handle
        self.MiPSHeater_button = MiPSHeater_button
    
    def run(self):
         # check and block the communication with MiPS while reading
         global MiPS_communication
         # check the communication is free, if not, wait till it is freed
         while MiPS_communication > 0:
             time.sleep(0.05)
         # block the communication
         MiPS_communication = 1
         # read current state of heater
         print('Not yet functional')
#         if str(self.MiPS_handle.getMiPSHeater) == 'OFF':
#             self.MiPS_handle.setMiPSHeater(self.MiPS_handle.devices['Magnet'], 'ON')
#             # free the communication
#             MiPS_communication = 0
#             self.MiPSHeater_button.setText('Heater ... on')
#             time.sleep(30)
#             self.MiPSHeater_button.setText('Heater ON')
#         elif str(self.MiPS_handle.getMiPSHeater) == 'ON':
#             self.MiPS_handle.setMiPSHeater(self.MiPS_handle.devices['Magnet'], 'OFF')
#             # free the communication
#             MiPS_communication = 0
#             self.MiPSHeater_button.setText('Heater ON')
         MiPS_communication = 0
             
class sub_read_current_field(QRunnable):
    def __init__(self, MiPS_handle, cFieldT):
        super().__init__()
        self.MiPS_handle = MiPS_handle
        self.cFieldT = cFieldT
        self.cFieldT_value = self.cFieldT.value()
        
    def run(self):
         # check and block the communication with MiPS while reading
         global MiPS_communication
         # check the communication is free, if not, wait till it is freed
         while MiPS_communication > 0:
             time.sleep(0.05)
         # block the communication
         MiPS_communication = 1
         # Read current field strength
         """ Read value here """
         self.cFieldT_value = self.MiPS_handle.getField('Magnet')
         self.cFieldT.display(float(self.cFieldT_value))
         # free the communication
         MiPS_communication = 0
         print('Read field done!')
    
class sub_set_field(QRunnable):
    def __init__(self, MiPS_handle, setMagnetSetPoint_box):
        super().__init__()
        self.MiPS_handle = MiPS_handle
        self.setMagnetSetPoint_box = setMagnetSetPoint_box
        
    def run(self):
         # check and block the communication with MiPS while reading
         global MiPS_communication
         # check the communication is free, if not, wait till it is freed
         while MiPS_communication > 0:
             time.sleep(0.05)
         # block the communication
         MiPS_communication = 1
         # Set field
         self.MiPS_handle.setField('Magnet', self.setMagnetSetPoint_box.value())
         # free the communication
         MiPS_communication = 0
         print('Set field done!')

class sub_set_ramp(QRunnable):
    def __init__(self, MiPS_handle, setMagnetRamp_box):
        super().__init__()
        self.MiPS_handle = MiPS_handle
        self.setMagnetRamp_box = setMagnetRamp_box
        
    def run(self):
         # check and block the communication with MiPS while reading
         global MiPS_communication
         # check the communication is free, if not, wait till it is freed
         while MiPS_communication > 0:
             time.sleep(0.05)
         # block the communication
         MiPS_communication = 1
         # Set field
         self.MiPS_handle.setFieldRate('Magnet', self.setMagnetRamp_box.value())
         # free the communication
         MiPS_communication = 0
         print('Set ramp done!')
         
class sub_MiPS_GoToSetPoint(QRunnable):
    def __init__(self, MiPS_handle):
        super().__init__()
        self.MiPS_handle = MiPS_handle
        
    def run(self):
         # check and block the communication with MiPS while reading
         global MiPS_communication
         # check the communication is free, if not, wait till it is freed
         while MiPS_communication > 0:
             time.sleep(0.05)
         # block the communication
         MiPS_communication = 1
         # Set field
         self.MiPS_handle.setAction('Magnet', 'RTOS')
         # free the communication
         MiPS_communication = 0
         print('Go to set point...')
         
class sub_MiPS_ToZero(QRunnable):
    def __init__(self, MiPS_handle):
        super().__init__()
        self.MiPS_handle = MiPS_handle
        
    def run(self):
         # check and block the communication with MiPS while reading
         global MiPS_communication
         # check the communication is free, if not, wait till it is freed
         while MiPS_communication > 0:
             time.sleep(0.05)
         # block the communication
         MiPS_communication = 1
         # Set field
         self.MiPS_handle.setAction('Magnet', 'RTOZ')
         # free the communication
         MiPS_communication = 0
         print('Go to zero...')
         
class sub_MagnetHold(QRunnable):
    def __init__(self, MiPS_handle):
        super().__init__()
        self.MiPS_handle = MiPS_handle
        
    def run(self):
         # check and block the communication with MiPS while reading
         global MiPS_communication
         # check the communication is free, if not, wait till it is freed
         while MiPS_communication > 0:
             time.sleep(0.05)
         # block the communication
         MiPS_communication = 1
         # Set field
         self.MiPS_handle.setAction('Magnet', 'HOLD')
         # free the communication
         MiPS_communication = 0
         print('Set on Hold')
# END: Subclasses for parallel processing


