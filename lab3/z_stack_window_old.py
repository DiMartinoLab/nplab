"""
Dawn - dmk50
October 2022 - aiming to have new guit with control of both ocean optics and the stage - to enable df-stack
Trung - xtn20
November 2022 - SmarAct stage control with MCS2 controller
"""
#nb, next up - look at spectrometer displayui

from PyQt5 import QtWidgets, uic
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

#import curve to be used when using reflective metal as reference material rather than white scatterer
ratio_curve_points = np.genfromtxt(fname='ratio_curve.txt')



class z_stack_window_object(Experiment, QtWidgets.QWidget, UiTools):
    #inherit QWidget methods
    def __init__(self,activeDatafile, OOSpect_instance, Stage_instance, stage_enabled,  ui_file = os.path.join(os.path.dirname(__file__),'z_stack_window.ui'),   parent=None):
        
        # set objects for spectrometer and stage
        self.OOspectrometer = OOSpect_instance
        self.stage_enabled = stage_enabled
        if self.stage_enabled == 'Cryostat stage':
            self.Smaract_stage = Stage_instance
        elif self.stage_enabled == 'Olympus stage':
            self.SMC_stage = Stage_instance

            
        # make this window available for the main program
        super(z_stack_window_object, self).__init__() 
        uic.loadUi(ui_file, self)
        
        
        self.z_stack_df_group = activeDatafile.create_group(name='stack group', auto_increment=True)
        #dset = self.z_stack_df_group.create_dataset(name, *args, **kwargs)
        self.z_stack_df_group.attrs.create('test attribute', str(101))
        #activeDatagroup = activeDatafile.create_group('z_stack_data_group', auto_increment=True,)
        #activeDatagroup.attrs.create('stage_enabled ', self.stage_enabled )        
        #df.create_group(name, auto_increment=True, *args, **kwargs)
        
        self.ratio_curve_values = [item[1] for item in ratio_curve_points]
        """connect spectrometer temperature control and readout buttons"""
        self.set_tec_temperature_pushButton.clicked.connect(self.gui_set_tec_temperature)
        self.read_tec_temperature_pushButton.clicked.connect(self.gui_read_tec_tempeature)
        initial_temperature = np.round(self.OOspectrometer.get_tec_temperature(), decimals = 1)
        self.tec_temperature_lcdNumber.display(float(initial_temperature))
        self.set_tec_temperature_LineEdit.setText(str(initial_temperature))
        """connect spectrometer integration time control and readout buttons - nb, may need to check second - millisecond conversion here"""    
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
        self.averaging_enabled=False
        
        """setup display for spectrum plot """
        self.plotbox = QtWidgets.QGroupBox()
        self.plotbox.setLayout(QtWidgets.QGridLayout())
        self.plotlayout = self.plotbox.layout()     
        self.spec_plot=pg.PlotWidget(labels = {'bottom':'Wavelength (nm)'})
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
            channel = 0
            smaract_package.SetProperty_i32(self.Smaract_stage, channel, smaract_package.Property.MOVE_MODE, self.smaract_movemode)
            smaract_package.SetProperty_i32(self.Smaract_stage, channel, smaract_package.Property.STEP_FREQUENCY, 1000)
            smaract_package.SetProperty_i32(self.Smaract_stage, channel, smaract_package.Property.STEP_AMPLITUDE, 65535)
            channel = 1
            smaract_package.SetProperty_i32(self.Smaract_stage, channel, smaract_package.Property.MOVE_MODE, self.smaract_movemode)
            smaract_package.SetProperty_i32(self.Smaract_stage, channel, smaract_package.Property.STEP_FREQUENCY, 1000)
            smaract_package.SetProperty_i32(self.Smaract_stage, channel, smaract_package.Property.STEP_AMPLITUDE, 65535)
            channel = 2
            smaract_package.SetProperty_i32(self.Smaract_stage, channel, smaract_package.Property.MOVE_MODE, self.smaract_movemode)
            smaract_package.SetProperty_i32(self.Smaract_stage, channel, smaract_package.Property.STEP_FREQUENCY, 1000)
            smaract_package.SetProperty_i32(self.Smaract_stage, channel, smaract_package.Property.STEP_AMPLITUDE, 65535)
            
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
        self.OOspectrometer.reference = self.OOspectrometer.read_spectrum() 
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
        self.current_spec = self.OOspectrometer.read_spectrum()
        self.current_spec = self.process_spec(self.current_spec)
        
    def take_spec_and_update_display(self):
        self.current_spec = self.OOspectrometer.read_spectrum()
        self.current_spec = self.process_spec(self.current_spec)
        self.update_display(self.current_spec)
        
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
            
    def z_stack_method(self):
        self.stack_plot.clear()
    
        self.stack_step_size = int(self.zstep_box_2.value())
        self.stack_num_steps = int(self.num_steps_box.value())
        # use this line to move away from sample:
        step_num=0
        self.stack_spectra = []
        while step_num<self.stack_num_steps:
            self.take_spec()
            self.smaract_move(1, self.stack_step_size*(-1))
            self.stack_spectra = self.stack_spectra + [self.current_spec]
            self.stack_plot.plot(self.wl_in_nm, self.current_spec)
            step_num = step_num+1
        self.make_stack_plot()
            
    def make_stack_plot(self):
        spec_nums=np.arange(0,self.stack_num_steps,1)
        spectra_T=np.transpose(self.stack_spectra)
        im=plt.pcolormesh(spec_nums, self.wl_in_nm, spectra_T, shading='auto', norm=colors.SymLogNorm(base=np.e, linthresh=0.5, linscale=0.05, vmin=-0.5, vmax=0.5))
        plt.clim(0,0.1)
        plt.show()
        
        
    # Trung Nov. 2022
    # added parts for the SmarAct stage control
    # Further infos and code: C:\SmarAct\MCS2\SDK\Python\examples
    # self.Smaract_stage = d_handle
  
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
        position = smaract_package.GetProperty_i64(self.Smaract_stage, 1, smaract_package.Property.POSITION)
        print("MCS2 position: {} pm.".format(position))
    # added parts for the SMC stage control
    
    
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


"""

    def save_spectrum(self, spectrum=None, attrs={}, new_deque = False):
        #Save a spectrum to the current datafile, creating if necessary.
        
        #If no spectrum is passed in, a new spectrum is taken.  The convention
        #is to save raw spectra only, along with reference/background to allow
        #later processing.
        
        #The attrs dictionary allows extra metadata to be saved in the HDF5 file
        if self.averaging_enabled == True:
            spectrum = self.read_averaged_spectrum(new_deque = new_deque)
        else:
            spectrum = self.read_spectrum() if spectrum is None else spectrum
        metadata = self.metadata
        metadata.update(attrs) #allow extra metadata to be passed in
        self.create_dataset(self.filename, data=spectrum, attrs=metadata) 
        #save data in the default place (see nplab.instrument.Instrument)
        

        
  
"""

