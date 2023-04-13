"""
Dawn - dmk50
October 2022 - aiming to have new gui with control of both ocean optics and the stage - to enable df-stack
April 2023   - adding option to use SMC100 stage rather than the cryostat smaract stage 

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
        
        # set objects for spectrometer 
        self.OOspectrometer = OOSpect_instance
        self.OOspectrometer.set_integration_time(50)
        self.set_default_int_times( default_t_in_ms=50)
        global OO_spectrometer_in_use
        OO_spectrometer_in_use = False
        
        #check wch stage is connected
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
        self.time_series_group = self.StackGUI_df_group.create_group(name='spectrum time series ', auto_increment=True)
        self.ratio_curve_values = [item[1] for item in ratio_curve_points]
        self.OO_Spectra_group.attrs.create('Ratio_curve', self.ratio_curve_values)

        
        """connect control and readout buttons for spectrometer temperature & spectrometer integration time """
        self.set_tec_temperature_pushButton.clicked.connect(self.gui_set_tec_temperature)
        self.read_tec_temperature_pushButton.clicked.connect(self.gui_read_tec_tempeature)
        initial_temperature = np.round(self.OOspectrometer.get_tec_temperature(), decimals = 1)
        self.tec_temperature_lcdNumber.display(float(initial_temperature))
        self.set_tec_temperature_LineEdit.setText(str(initial_temperature))
        self.set_integration_time_push_button.clicked.connect(self.gui_set_integration_time)
        
        """setup spectrometer commands """
        self.wl_in_nm = self.OOspectrometer.get_wavelengths()
        
        self.referencing_type = 'Gold'
        self.white_ref_checkbox.setChecked(False)
        self.gold_ref_checkbox.setChecked(True)
        self.white_ref_checkbox.stateChanged.connect(self.reference_choice_toggle)
        self.gold_ref_checkbox.stateChanged.connect(self.reference_choice_toggle)
        
        self.read_bkg_for_ref_button.clicked.connect(self.read_bkg_for_ref_method)
        self.clear_bkg_for_ref_button.clicked.connect(self.clear_bkg_for_ref_method)
        self.read_bkg_button.clicked.connect(self.read_background_method)
        self.clear_bkg_button.clicked.connect(self.clear_background_method)
        self.read_ref_button.clicked.connect(self.read_ref_method)
        self.clear_ref_button.clicked.connect(self.clear_ref_method)
        
        self.take_spectrum_button.clicked.connect(self.take_spec_and_update_display)
        self.save_button.clicked.connect(self.save_single_DF_spectrum)
        
        self.time_series_pushButton.clicked.connect(self.take_t_series)
    
        self.OO_disconnect_button.clicked.connect(self.OO_disconnect)
        
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
        
        """setup stage controls and display - using olympus OR smaract"""
        self.disconnect_stage_button.clicked.connect(self.disconnect_stage)
        # dmk50 April 2023
        self.stage_mleft_button.clicked.connect(self.SMC100_move_left)
        if stage_enabled == 'Olympus stage':
            #this is used to prevent sending too many commands to stage at once
            global stage_in_use
            stage_in_use = False
            
            #make and connect button to initialise the stage
            self.stage_feature_button.setText('Init SMC100 ')
            self.stage_feature_button.clicked.connect(self.initialise_SMC100)
            
            # link buttons for stage movement to functions
            self.stage_mleft_button.clicked.connect(self.SMC100_move_left)
            self.stage_mright_button.clicked.connect(self.SMC100_move_right)
            self.stage_mup_button.clicked.connect(self.SMC100_move_up)
            self.stage_mdown_button.clicked.connect(self.SMC100_move_down)
            self.stage_mtoward_button.setText('Towards lens')
            self.stage_mtoward_button.clicked.connect(self.SMC100_move_towardslens)
            self.stage_maway_button.setText('Away from lens')
            self.stage_maway_button.clicked.connect(self.SMC100_move_awaylens)
            self.stage_mmid_button.clicked.connect(self.SMC100_move_mid)
            #setup button to print stage position 
            self.get_stage_position_button.clicked.connect(self.print_stage_positions)
            
        
        # Trung Nov. 2022
        # Further infos and code: C:\SmarAct\MCS2\SDK\Python\examples
        if stage_enabled == 'Cryostat stage':
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
            self.stage_feature_button.setText('Camera view')
            
            self.smaract_channel_up = 0
            self.smaract_direction_up = -1
            
            self.smaract_channel_down = 0
            self.smaract_direction_down = 1
            
            self.smaract_channel_left = 2
            self.smaract_direction_left = -1
            
            self.smaract_channel_right = 2
            self.smaract_direction_right = 1


            # link buttons for smaract movement to functions
            self.stage_mleft_button.clicked.connect(self.smaract_moveleft)
            self.stage_mright_button.clicked.connect(self.smaract_moveright)
            self.stage_mup_button.clicked.connect(self.smaract_moveup)
            self.stage_mdown_button.clicked.connect(self.smaract_movedown)
            
            self.stage_mtoward_button.clicked.connect(self.smaract_movetoward)
            self.stage_maway_button.clicked.connect(self.smaract_moveaway)
        
            # Movement of stage on cryostat may differ from movement on camera view
            # Use stage feature button to change view - and display view name
            self.stage_feature_button.clicked.connect(self.current_view)
        
        
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
        
        
    """START: methods for use with OO spectrometer"""
    
    def set_default_int_times(self, default_t_in_ms):
        self.OOspectrometer.integration_time  = default_t_in_ms
        self.OOspectrometer.background_int_ref= default_t_in_ms
        self.OOspectrometer.background_int    = default_t_in_ms
        self.OOspectrometer.reference_int     = default_t_in_ms

    def read_bkg_for_ref_method(self):
        #read spec and store as background for reference. chenge check-box state to match.
        self.OOspectrometer.background_ref = self.OOspectrometer.read_spectrum()
        self.OOspectrometer.background_int_ref = self.OOspectrometer.integration_time
        self.OOspectrometer.stored_backgrounds[self.OOspectrometer.reference_ID] = {'background_ref' : self.OOspectrometer.background_ref,'background_for_ref_int_t': self.OOspectrometer.background_int_ref}
        self.OOspectrometer.update_config('background_ref', self.OOspectrometer.background_ref)
        self.OOspectrometer.update_config('background_for_ref_int_t', self.OOspectrometer.background_int_ref)
        print(self.OOspectrometer.read_background_ref)
        self.background_subtracted_ref.setCheckState(QtCore.Qt.Checked)
            
    def clear_bkg_for_ref_method(self):
        self.OOspectrometer.background_ref = None
        self.background_subtracted_ref.setCheckState(QtCore.Qt.Unchecked)
        
    def fetch_bkg_for_ref_method(self, a_spectrum):
        # if no bkg for ref has been collected, set to array of zeros of appropriate length
        if self.background_subtracted_ref.isChecked()==True:
            pass
        elif self.background_subtracted_ref.isChecked()==False:
            self.OOspectrometer.background_ref = np.zeros(len(a_spectrum))
            
    def read_ref_method(self):
        self.OOspectrometer.reference = self.OOspectrometer.read_spectrum() 
        self.OOspectrometer.reference_int = self.OOspectrometer.integration_time
        self.OOspectrometer.update_config('reference', self.OOspectrometer.reference)
        self.OOspectrometer.update_config('reference_int',self.OOspectrometer.reference_int ) 
        self.OOspectrometer.stored_references[self.OOspectrometer.reference_ID] = {'reference' : self.OOspectrometer.reference,'reference_int' : self.OOspectrometer.reference_int}
        self.referenced.setCheckState(QtCore.Qt.Checked)

    def clear_ref_method(self):
        self.OOspectrometer.reference = None
        self.referenced.setCheckState(QtCore.Qt.Unchecked)
        
    def fetch_ref_method(self, a_spectrum):
        # if no  ref has been collected, set to array of ones of appropriate length
        if self.referenced.isChecked()==True:
            pass
        elif self.referenced.isChecked()==False:
            self.OOspectrometer.reference = np.ones(len(a_spectrum))

    def read_background_method(self):
        """Acquire a new spectrum and use it as a background measurement.
        This background should be less than 50% of the spectrometer saturation"""
        self.OOspectrometer.background = self.OOspectrometer.read_spectrum()
        self.OOspectrometer.background_int = self.OOspectrometer.integration_time
        self.OOspectrometer.stored_backgrounds[self.OOspectrometer.reference_ID] = {'background' : self.OOspectrometer.background,
                                                     'background_int': self.OOspectrometer.background_int}
        self.OOspectrometer.update_config('background', self.OOspectrometer.background)
        self.OOspectrometer.update_config('background_int', self.OOspectrometer.background_int)
        self.background_subtracted.setCheckState(QtCore.Qt.Checked)
        
    def clear_background_method(self):
        self.OOspectrometer.background = None
        #self.OOspectrometer.background_int = None
        self.background_subtracted.setCheckState(QtCore.Qt.Unchecked)
    
    def fetch_bkg_method(self, a_spectrum):
        # if no  bkg has been collected, set to array of ones of appropriate length
        if self.background_subtracted.isChecked()==True:
            pass
        elif self.background_subtracted.isChecked()==False:
            self.OOspectrometer.background  = np.zeros(len(a_spectrum))
    
    def take_spec(self):

        #update int t
       
        # read the current field strength and sample temperature
        if self.MiPS_handle != None:
            current_field = sub_read_current_field(self.MiPS_handle, self.cFieldT)
            self.pool.start(current_field)
        if self.MiTC_handle  != None:
            current_temp = sub_read_current_sample_temp(self.MiTC_handle, self.cSampleTemp)
            self.pool.start(current_temp)
        
        #hold cryostat stage steady during integration of signal by the spectrometer
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
            
        elif self.stage_enabled == 'Olympus stage':
            self.current_z_position = self.SMC_stage.get_position(3)

        
        # Take the spectrum in the background
        global OO_spectrometer_in_use
        OO_spectrometer_in_use=True
        take_spec_task = sub_take_spec(self.OOspectrometer)
        self.pool.start(take_spec_task)
        
        # Wait till spectrum is taken
        
        while OO_spectrometer_in_use ==True:
            pass
        print(OO_spectrometer_in_use )
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
            
        #storing raw and processed spec in seperate variables
        self.current_processed_spec = self.process_spec(self.current_raw_spec)
        
    def take_spec_and_update_display(self):
        # take_spec method allows parallel processing
        self.take_spec()
        #this line in twice to get rid of weird timing glitch...temporary patch...
        self.take_spec()
        global OO_spectrometer_in_use 
        while OO_spectrometer_in_use == True:
            pass
        if OO_spectrometer_in_use == False:
            self.update_display(self.current_processed_spec)
        
    
        
    def process_spec(self, a_spectrum):
        # Fetch bkg for ref, ref, and bkg
        # If any of above not stored - return array of zeros or ones as appropriate to stand as placeholder in calculation
        self.fetch_bkg_for_ref_method(a_spectrum)
        self.fetch_bkg_method(a_spectrum)
        
        self.fetch_ref_method(a_spectrum)
        
        #scale bkg/ref spectra according to ratio of integration times with data spectrum being processed
        scaled_bkg_for_ref = self.OOspectrometer.background_ref * (self.OOspectrometer.integration_time/self.OOspectrometer.background_int_ref)
        scaled_bkg = self.OOspectrometer.background * (self.OOspectrometer.integration_time/self.OOspectrometer.background_int)
        scaled_ref = self.OOspectrometer.reference * (self.OOspectrometer.integration_time/self.OOspectrometer.reference_int)
        
        #check whether white or metal reference is used - use/don't use ratio curve accordingly
        #nb - only use ratio curve if a reference spectrum has been measured
        if self.referencing_type == 'Gold':
            if self.referenced.isChecked()==False:
                scaled_ref_bkg_subtract = (scaled_ref - scaled_bkg_for_ref )
            elif self.referenced.isChecked()==True:
                scaled_ref_bkg_subtract = (scaled_ref - scaled_bkg_for_ref )* self.ratio_curve_values 
        elif self.referencing_type == 'White':
            scaled_ref_bkg_subtract = (scaled_ref - scaled_bkg_for_ref )
            
        # accounting for case where bkg_for_ref measured but ref not: 
        if self.referenced.isChecked()==False and self.background_subtracted_ref.isChecked()==True:
            scaled_ref_bkg_subtract = np.ones(len(a_spectrum))
            
        #calculate processed spectrum
        processed_spectrum = (a_spectrum - scaled_bkg)/(scaled_ref_bkg_subtract)
        return processed_spectrum 
        

    def update_display(self, spectrum):
        #replace any Nan in input spectrum with a 0
        spectrum = np.array([0 if np.isnan(i) else i for i in spectrum])
        #empty the display each time so that only most recent plot is shown
        self.spec_plot.clear()
        self.spec_plot.plot(self.wl_in_nm, spectrum)

    """connect spectrometer control and readout buttons"""
    def gui_set_tec_temperature(self):
        self.OOspectrometer.tec_temperature = float(self.set_tec_temperature_LineEdit.text().strip())
        
    def gui_read_tec_tempeature(self):
        self.tec_temperature_lcdNumber.display(float(self.OOspectrometer.tec_temperature)) 
        
    def gui_set_integration_time(self):
        self.OOspectrometer.integration_time = float(self.set_integration_time_LineEdit.text().strip())
        self.OOspectrometer.set_integration_time(self.OOspectrometer.integration_time)
        print('Intergration time set to: '+str(self.OOspectrometer.integration_time))
        
    def reference_choice_toggle(self, state):
        sender = self.sender()
        if sender is self.gold_ref_checkbox and state == QtCore.Qt.Checked:
            print('Using Gold Ref')
            self.white_ref_checkbox.setChecked(False)
            self.referencing_type = 'Gold'
        elif sender is self.white_ref_checkbox and state == QtCore.Qt.Checked:
            print('Using White Ref')
            self.gold_ref_checkbox.setChecked(False)
            self.referencing_type = 'White'
            
    def read_spec_description(self):
            self.current_spec_description= str(self.description_LineEdit.text().strip())
            
    def save_single_DF_spectrum(self):
        #save a raw spectrum = self.current_raw_spec
        """nb - want to save both raw and processed!"""
        activeSingleDFDataset = self.OO_Spectra_group.create_dataset('singleDFSpectrum_%d', data = self.current_raw_spec)

        self.read_spec_description()
        activeSingleDFDataset.attrs.create("Description", self.current_spec_description)
        activeSingleDFDataset.attrs.create("Referencing Type", self.referencing_type)
        activeSingleDFDataset.attrs.create("Processed Spectrum", self.current_processed_spec)
        activeSingleDFDataset.attrs.create("Wavelengths", self.wl_in_nm)
        activeSingleDFDataset.attrs.create('Integration_Time', self.OOspectrometer.integration_time)
    
        activeSingleDFDataset.attrs.create('bkg_for_ref_int_t', self.OOspectrometer.background_int_ref)
        activeSingleDFDataset.attrs.create('bkg_int_t', self.OOspectrometer.background_int)
        activeSingleDFDataset.attrs.create('ref_int_t', self.OOspectrometer.reference_int)
        
        activeSingleDFDataset.attrs.create('bkg_for_ref', self.OOspectrometer.background_ref )
        activeSingleDFDataset.attrs.create('background', self.OOspectrometer.background)
        activeSingleDFDataset.attrs.create('reference', self.OOspectrometer.reference)
        
        activeSingleDFDataset.attrs.create('solenoid temperature', self.current_temp)
        activeSingleDFDataset.attrs.create('solenoid field', self.current_field)
        
        activeSingleDFDataset.attrs.create('stage z position', self.current_z_position )
        print(' ')
        print('A single DF spectrum has been saved.' )
        print('Solenoid params =  '+str(self.current_field) + ' T and ' + str(self.current_temp) + ' K.')        
        
        
    def take_t_series(self):
        # num_spectra_spinBox = box containing num spectra to be taken in t-series
        # time_delay_SpinBox =  box containing time delay to be between spectra
        # time_series_name_lineEdit = text box containing name to save t-series under
        t_series_num_spec = int(self.num_spectra_spinBox.value())
        t_series_delay    = float(self.time_delay_SpinBox.value())/1000
        t_series_name     = str(self.time_series_name_lineEdit.text().strip())
        self.raw_t_series_spectra       = []
        self.processed_t_series_spectra =[]
        self.t_series_times = []
        
        spectrum_counter = 0
        time_counter = 0
        print('  BEGIN TIME SERIES')
        print('t_series_delay = '+str(t_series_delay*1000)+' ms.')
        print(' self.OOspectrometer.integration_time = '+str(self.OOspectrometer.integration_time)+' ms.')
        while spectrum_counter < t_series_num_spec:
            self.take_spec()
            self.t_series_times             = self.t_series_times             + [time_counter ]
            self.raw_t_series_spectra       = self.raw_t_series_spectra       + [self.current_raw_spec]
            self.processed_t_series_spectra = self.processed_t_series_spectra + [self.current_processed_spec]
            time.sleep(t_series_delay)
            time_counter = time_counter + self.OOspectrometer.integration_time + (t_series_delay*1000)
            print('spectrum_counter   = '+str(spectrum_counter ))
            spectrum_counter  = spectrum_counter + 1

        activeTSeriesDataset = self.time_series_group.create_dataset('SpectrumTimeSeries_%d', data = self.raw_t_series_spectra)
        
        activeTSeriesDataset.attrs.create('t-series num spectra', t_series_num_spec)
        activeTSeriesDataset.attrs.create('t-series spec end-to-start delay', t_series_delay)
        activeTSeriesDataset.attrs.create("Wavelengths", self.wl_in_nm)
        activeTSeriesDataset.attrs.create('Integration_Time', self.OOspectrometer.integration_time)         
        activeTSeriesDataset.attrs.create("z-position", self.current_z_position)
        activeTSeriesDataset.attrs.create("T series name", t_series_name)
        #activeTSeriesDataset.attrs.create("Processed Spectra", self.processed_t_series_spectra )
        
        activeTSeriesDataset.attrs.create('bkg_for_ref_int_t', self.OOspectrometer.background_int_ref)
        activeTSeriesDataset.attrs.create('bkg_int_t', self.OOspectrometer.background_int)
        activeTSeriesDataset.attrs.create('ref_int_t', self.OOspectrometer.reference_int)
        activeTSeriesDataset.attrs.create('bkg_for_ref', self.OOspectrometer.background_ref )
        activeTSeriesDataset.attrs.create('background', self.OOspectrometer.background)
        activeTSeriesDataset.attrs.create('reference', self.OOspectrometer.reference)
        
        activeTSeriesDataset.attrs.create('solenoid temperature', self.current_temp)
        activeTSeriesDataset.attrs.create('solenoid field', self.current_field)
        print('  END TIME SERIES')
        
    """END: methods for use with OO spectrometer"""
            
    """START: methods for z-stacks using OO spectrometer"""
    def z_stack_method(self):
        self.stack_plot.clear()
        self.stack_step_size = int(self.zstep_box_2.value())
        self.stack_num_steps = int(self.num_steps_box.value())
        
        global z_stack_running
        z_stack_running = 1
        

        if self.stage_enabled == 'Cryostat stage':
            initial_z = smaract_package.GetProperty_i64(self.Smaract_stage, 1, smaract_package.Property.POSITION)
        elif self.stage_enabled == 'Olympus stage':
            initial_z =self.SMC_stage.get_position(3)
                
        self.processed_stack_spectra = []
        self.raw_stack_spectra = []
        step_num=0
        while step_num < self.stack_num_steps:
            self.take_spec()
            global finished_OO
            while finished_OO < 1:
                pass
            if finished_OO==1:
                if self.stage_enabled == 'Cryostat stage':
                    self.smaract_move(1, self.stack_step_size*(-1))
                elif self.stage_enabled == 'Olympus stage':
                    print(' ')
                    print('z stack implementation with olympus stage not yet tested')
                    self.SMC100_move_zaxis(self.stack_step_size)
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
        
        if self.stage_enabled == 'Cryostat stage':
            final_z = smaract_package.GetProperty_i64(self.Smaract_stage, 1, smaract_package.Property.POSITION)
        elif self.stage_enabled == 'Olympus stage':
            final_z=self.SMC_stage.get_position(3)
            
        
        activeStackDataset.attrs.create("Final z position", final_z)
        self.read_spec_description()
        activeStackDataset.attrs.create("Description", self.current_spec_description)
        activeStackDataset.attrs.create("Processed Spectra", self.processed_stack_spectra )
        
        activeStackDataset.attrs.create('bkg_for_ref_int_t', self.OOspectrometer.background_int_ref)
        activeStackDataset.attrs.create('bkg_int_t', self.OOspectrometer.background_int)
        activeStackDataset.attrs.create('ref_int_t', self.OOspectrometer.reference_int)
        activeStackDataset.attrs.create('bkg_for_ref', self.OOspectrometer.background_ref )
        activeStackDataset.attrs.create('background', self.OOspectrometer.background)
        activeStackDataset.attrs.create('reference', self.OOspectrometer.reference)
        
        activeStackDataset.attrs.create('solenoid temperature', self.current_temp)
        activeStackDataset.attrs.create('solenoid field', self.current_field)
        self.make_stack_plot()

        #TRUNG TO DO TEST IPS
            
    def make_stack_plot(self):
        spec_nums = np.arange(0,self.stack_num_steps,1)
        spectra_T = np.transpose(self.processed_stack_spectra)
        plt.pcolormesh(spec_nums, self.wl_in_nm, spectra_T, shading = 'auto', norm = colors.SymLogNorm(base = np.e, linthresh = 0.5, linscale = 0.05, vmin = -0.5, vmax = 0.5))
        plt.clim(0,0.1)
        plt.show()
    """END: methods for df spectrum z-stacks"""
    
    
        
        
        
    # Trung Nov. 2022
    # added parts for the SmarAct stage control
    # Further infos and code: C:\SmarAct\MCS2\SDK\Python\examples
    # self.Smaract_stage = d_handle
  
    
    
    """START: methods for control of stages"""
    def disconnect_stage(self):
        if self.stage_enabled == 'Cryostat stage':
            smaract_package.Close(self.Smaract_stage)
            print('Smaract stage disconnected.')
        elif self.stage_enabled == 'Olympus stage':
            self.SMC_stage.__del__()
            print('Olynpus stage disconnected.')
            
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
        button_text = self.stage_feature_button.text()
        if button_text == 'Camera view':
            self.stage_feature_button.setText('Cryostat view')
            
            self.smaract_channel_left = 0
            self.smaract_direction_left = 1
            
            self.smaract_channel_right = 0
            self.smaract_direction_right = -1
            
            self.smaract_channel_up = 2
            self.smaract_direction_up = -1
            
            self.smaract_channel_down = 2
            self.smaract_direction_down = 1
        
        elif button_text == 'Cryostat view':
            self.stage_feature_button.setText('Camera view')
            
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

    def initialise_SMC100(self):
        do_initialise = sub_initialise_stage(self.SMC_stage)
        self.pool.start(do_initialise)
        
    def SMC100_move_left(self):
        axis = 1
        move_amount = float(self.xystep_Box.value())*100/1e6
        do_move = sub_SMC100_move(self.SMC_stage, str(axis), move_amount)
        self.pool.start(do_move)
    
    def SMC100_move_right(self):
        axis = 1
        goto_value = -float(self.xystep_Box.value())*100/1e6
        do_move = sub_SMC100_move(self.SMC_stage, str(axis), goto_value)
        self.pool.start(do_move)
    
    def SMC100_move_up(self):
        axis = 2
        goto_value = float(self.xystep_Box.value())*100/1e6
        do_move = sub_SMC100_move(self.SMC_stage, str(axis), goto_value)
        self.pool.start(do_move)
    
    def SMC100_move_down(self):
        axis = 2
        goto_value = -float(self.xystep_Box.value())*100/1e6
        do_move = sub_SMC100_move(self.SMC_stage, str(axis), goto_value)
        self.pool.start(do_move)

    def SMC100_move_towardslens(self):
        axis = 3
        goto_value = float(self.zstep_box.value())*50/1e6
        do_move = sub_SMC100_move(self.SMC_stage, str(axis), goto_value)
        self.pool.start(do_move)
    
    def SMC100_move_awaylens(self):
        axis = 3
        goto_value = -float(self.zstep_box.value())*50/1e6
        do_move = sub_SMC100_move(self.SMC_stage, str(axis), goto_value)
        self.pool.start(do_move)
        
    def SMC100_move_zaxis(self, step_size):
        axis = 3
        goto_value = -float(step_size)*50/1e6
        do_move = sub_SMC100_move(self.SMC_stage, str(axis), goto_value)
        self.pool.start(do_move)
        
    def SMC100_move_mid(self):
        do_move_mid = sub_SMC100_move_mid(self.SMC_stage)
        self.pool.start(do_move_mid)
        
    def print_stage_positions(self):
        #note - cryostat stage case is untested
        #note - need to add units to position values
        
        if self.stage_enabled == 'Cryostat stage':
            print('Smaract stage positions: current view = '+str(self.stage_feature_button.text()))
            z_position = smaract_package.GetProperty_i64(self.Smaract_stage, 1, smaract_package.Property.POSITION)
            left_right_position = smaract_package.GetProperty_i64(self.Smaract_stage, self.smaract_channel_left, smaract_package.Property.POSITION)
            up_down_position = smaract_package.GetProperty_i64(self.Smaract_stage, self.smaract_channel_up, smaract_package.Property.POSITION)
            print('     - - - - -     ')
            print('left/right axis = '+str(left_right_position)+' (moves in ?nm increments)')
            print('up/down axis = '+str(up_down_position)+' (moves in ?nm increments)')
            print('z-axis = '+str(z_position)+' (moves in ?nm increments)')
            print('     - - - - -     ')
            
        elif self.stage_enabled == 'Olympus stage':
            print('          ')
            print('Olympus stage positions:')
            print('     - - - - -     ')
            print('left/right axis = '+str(self.SMC_stage.get_position(1))+'mm (moves in 100nm increments)')
            print('up/down axis = '+str(self.SMC_stage.get_position(2))+'mm (moves in 100nm increments)')
            print('z-axis = '+str(self.SMC_stage.get_position(3))+'mm (moves in 50nm increments)')
            print('     - - - - -     ')
            

        
    """END: methods for control of stages"""
        
        
    """start: methods for control of cryostat field and temperatures"""
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
    """end: methods for control of cryostat field and temperatures"""
        
     
    def OO_disconnect(self):
        self.OOspectrometer._close()
        print('Ocean Optics successfully disconnected')
        

        
        

    def make_window(self):
        app = get_qt_app()
        self.show()
        app.exec_()
        return self













"""helper class for spectrometer"""
class sub_take_spec(QRunnable):
    #this class uses a seperate thread while taking a spectrum.
    def __init__(self, OOspectrometer):
        super().__init__()
        # this is local access to OOspectrometer. "self" here is the "sub_take_spec"
        self.OOspectrometer = OOspectrometer
    
    def run(self):
        self.current_spec = self.OOspectrometer.read_spectrum()
        global OO_spectrometer_in_use 
        OO_spectrometer_in_use = False
        print('Finished taking spec, using int t = '+str(self.OOspectrometer.integration_time))
    
"""start: helper classes for SMC100 stage"""
class sub_SMC100_move(QRunnable):
    def __init__(self, SMC100, axis, value):
        #axis is a string and value is absolute position 
        super().__init__()
        self.SMC100 = SMC100
        self.axis = axis
        self.value = value        
        global stage_in_use
        
    def run(self):
        global stage_in_use
        if stage_in_use ==True:
            print('Wait for stage release')
        while stage_in_use ==True:            
            time.sleep(0.5)
        print(' ')
        print('Moving SMC100')
        stage_in_use = True
        pos_stage = self.SMC100.get_position(self.axis)
        goto_pos = float(pos_stage[0]) + self.value
        self.SMC100.move(goto_pos, self.axis, relative=False)
        pos_stage = self.SMC100.get_position(self.axis)
        print('Position axis ' + str(self.axis) + ': ' + str(pos_stage) + ' mm')
        stage_in_use = False 
        

        
        
class sub_SMC100_move_mid(QRunnable):
    def __init__(self, SMC100):
        super().__init__()
        self.SMC100 = SMC100

    def run(self):
        global stage_in_use
        if stage_in_use ==True:
            print('Wait for stage release')
        while stage_in_use ==True:            
            time.sleep(0.5)
            print(' ')
        print('Moving to midpoint of stage range of motion.')
        stage_in_use = True
        
        self.SMC100.move(6, '1', relative=False)
        self.SMC100.move(6, '2', relative=False)
        self.SMC100.move(6, '3', relative=False)
        print('Stage all axes moved to position 6 mm')
        print(' ')
        stage_in_use = False 
        
        
class sub_initialise_stage(QRunnable):
    def __init__(self, SMC100):
        super().__init__()
        self.SMC100 = SMC100

    def run(self):
        self.SMC100.reset_and_configure()
        self.SMC100.home()
        print(' ')
        print('Stage initialised successfully')
        print(self.SMC100.get_position(1))
        print(self.SMC100.get_position(2))
        print(self.SMC100.get_position(3))
        print(' ')
"""end: helper classed for SMC100 stage"""
    
"""start: helper class for smaract stage"""    
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
        global OO_spectrometer_in_use
        while OO_spectrometer_in_use==True:
            cposition = smaract_package.GetProperty_i64(self.Smaract_stage, self.channel, smaract_package.Property.POSITION)
            if np.abs(cposition - self.position) > 500000: # = 500nm
                smaract_package.Move(self.Smaract_stage, self.channel, self.position, 0)
                print('Stage position corrected')
"""end: helper class for smaract stage"""  

"""start: helper classes for solenoid temperature controller"""  
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
"""end: helper classes for solenoid temperature controller"""               
         
         
"""start: helper classes for solenoid field controller"""  
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
"""end: helper classes for solenoid field controller"""  

