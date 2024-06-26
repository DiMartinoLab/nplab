"""
Dawn - dmk50
October 2022 - aiming to have new gui with control of both ocean optics and the stage - to enable df-stack
April 2023   - adding option to use SMC100 stage rather than the cryostat smaract stage 

Trung - xtn20
November 2022 - SmarAct stage control with MCS2 controller

Trung
April 2024 - merge z_stack_window and RamanSpectrometer to enable x-y scan of
Raman spectrum
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
from nplab.instrument.stage.SMC100_lib_zstack import SMC100 as SMC_package
from nplab.ui.ui_tools import UiTools
import ctypes
from ctypes import byref, c_int, c_ulong, c_double
import pyqtgraph as pg
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg
from matplotlib.figure import Figure
from PySide6.QtCharts import QChart, QChartView, QLineSeries, QDateTimeAxis, QValueAxis
from datetime import datetime
import time


#import curve to be used when using reflective metal as reference material rather than white scatterer
ratio_curve_points = np.genfromtxt(fname='ratio_curve.txt')

class xyz_stack_DF_Raman(Experiment, QtWidgets.QWidget, UiTools):
    #inherit QWidget methods
    def __init__(self, activeDatafile, \
                 OOSpect_instance, Stage_instance, MiTC_instance, MiPS_instance, \
                 stage_enabled, Kymera_handle, Andor_handle, threadCount, pool, \
                 ui_file = os.path.join(os.path.dirname(__file__),'xyz_stack_DF_Raman.ui'),   parent=None, \
                 pixel_number = 1024, pixel_width = 26, use_shifts = False, \
                 laser_wl = 782.5, white_shutter = None, settings_filepath = None, \
                 camera_index = None):
        
        # make this window available for the main program
        super(xyz_stack_DF_Raman, self).__init__() 
        uic.loadUi(ui_file, self)
        
        
        # set objects for spectrometer 
        self.OOspectrometer = OOSpect_instance
        self.OOspectrometer.set_integration_time(50)
        self.set_default_int_times(default_t_in_ms=50)
        self.OOspectrometer.OO_occupied = 0
        
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
        
        
        #name a group within the h5 file to save collected data into
        self.StackGUI_df_group = activeDatafile.create_group(name='Z-Stack_DataGroup', auto_increment=True)
        self.StackGUI_df_group.attrs.create('Version Date', str('Date of latest edit of Stack Gui = 21-11-2022'))
        #save spectra here, with addintional info such as stage location and cryostat params also saved.
        self.OO_Spectra_group= self.StackGUI_df_group.create_group(name='OO Spectra', auto_increment=True)
        self.Stacked_Spectra_group = self.StackGUI_df_group.create_group(name='Stacked Spectra', auto_increment=True)
        self.time_series_group = self.StackGUI_df_group.create_group(name='spectrum time series', auto_increment=True)
        self.ratio_curve_values = [item[1] for item in ratio_curve_points]
        self.OO_Spectra_group.attrs.create('Ratio_curve', self.ratio_curve_values)

        
        """connect control and readout buttons for spectrometer temperature & spectrometer integration time """
        self.read_tec_temperature_pushButton.clicked.connect(self.gui_read_tec_tempeature)

        self.set_integration_time_push_button.clicked.connect(self.gui_set_integration_time)
        self.OOspectrometer.tec_temperature = float(-20)
        
        """setup spectrometer commands """
        self.wl_in_nm = self.OOspectrometer.get_wavelengths()
        
        self.referencing_type = 'Gold'
        self.white_ref_checkbox.setChecked(False)
        self.gold_ref_checkbox.setChecked(True)
        self.white_ref_checkbox.stateChanged.connect(self.reference_choice_toggle)
        self.gold_ref_checkbox.stateChanged.connect(self.reference_choice_toggle)
        
        self.read_bkg_for_ref_button.clicked.connect(self.read_bkg_for_ref_method)
        self.clear_bkg_for_ref_button.clicked.connect(self.clear_bkg_for_ref_method)
        self.clear_bkg_for_ref_method()
        self.read_bkg_button.clicked.connect(self.read_background_method)
        self.clear_bkg_button.clicked.connect(self.clear_background_method)
        self.clear_background_method()
        self.read_ref_button.clicked.connect(self.read_ref_method)
        self.clear_ref_button.clicked.connect(self.clear_ref_method)
        self.clear_ref_method()
        
        self.take_spectrum_button.clicked.connect(self.take_spec_and_update_display)
        self.save_button.clicked.connect(self.save_single_DF_spectrum)
        
        self.RestoreDFsetting_button.clicked.connect(self.RestoreDFsetting)
        self.time_series_pushButton.clicked.connect(self.t_series_method)
        self.SetfigRange_button.clicked.connect(self.SetfigRange)
        self.OO_disconnect_button.clicked.connect(self.OO_disconnect)
        

        """setup z-stack controls and display """
        self.take_z_stack_button.clicked.connect(self.z_stack_method)

        """setup stage controls and display - using olympus OR smaract"""
        self.disconnect_stage_button.clicked.connect(self.disconnect_stage)
        
        # dmk50 April 2023
        
        #buttons and display for z setpoint
        self.Set_New_Coord_origin_button.clicked.connect(self.Set_New_Origin)
        self.Go_To_Stage_Origin_button.clicked.connect(self.Go_To_Stage_Origin)
        self.origin_stored_checkbox.setChecked(False)
            
        if stage_enabled == 'Olympus stage':
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
            
            #setup button to print stage position 
            self.get_stage_position_button.clicked.connect(self.print_stage_positions)
        
        
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
#        self.threadCount = QThreadPool.globalInstance().maxThreadCount()
#        self.pool = QThreadPool.globalInstance()
        self.threadCount = threadCount
        self.pool = pool
        
        """ Add part for Raman Trung 11.04.2024"""
        
        """ Initialise the spectrometer """
        self.kymera = Kymera_handle
        self.kymera.pixel_number = pixel_number
        self.kymera.pixel_width = pixel_width
        self.use_shifts = use_shifts
        self.white_shutter = white_shutter
        self.pixel_number = pixel_number
        
        # set grating and central wavelength, uncomment while use with main GUI
        self.kymera.SetGrating(grating_num=1)
        self.kymera.SetWavelength(800)
#        self.metadata_property_names += ('slit_width', 'wavelengths')
        self.ImageFlip = 0
        global wavelength_Raman
        wavelength_Raman = self.get_x_axis()
        self.kymera.wavelength = wavelength_Raman
        
        """ Initialise the camera """
        self.Andor = Andor_handle
        self.CurImage = None
        self.background = None
        self.backgrounded = False
        self.keep_shutter_open = False
        if settings_filepath is not None:
            self.load_params_from_file(settings_filepath)
        self.isAborted = False
        self.Andor.Andor_occupied = 0
        self.set_exposureTime_Raman()
        # Put the wavelength into Andor for SMU script and also for
        # sub-functions in this script since kymera is not frequently called in
        # the sub-functions
        self.Andor.laser_wl = laser_wl

        """ Buttons"""
        self.Centralwavelength_button_Raman.clicked.connect(self.set_centralwavelength_Raman)
        self.Exposure_button_Raman.clicked.connect(self.set_exposureTime_Raman)
        self.Capture_button_Raman.clicked.connect(self.CaptureSpectrum_Raman)
        self.Readmode_button_Raman.clicked.connect(self.Switch_read_mode_Raman)
        self.SetROI_button_Raman.clicked.connect(self.setROI_Raman)
        self.Abort_button_Raman.clicked.connect(self.Abort_Raman)
        self.Save_button_Raman.clicked.connect(self.Save_Raman)
        self.Repeat_button_Raman.clicked.connect(self.Capture_and_save_Raman)
        self.Live_button_Raman.clicked.connect(self.Live_mode_Raman)
        self.Andor_occupied_reset_button_Raman.clicked.connect(self.Andor_occupied_reset)
        
        self.Switch_scan_direction_button_Raman.clicked.connect(self.Raman_scan_direction)
        self.Capture_xy_scan_button_Raman.clicked.connect(self.Capture_xy_scan_Raman)
        global Abort_live_mode
        Abort_live_mode = 0
        

        self.Andor_temp_button_Raman.clicked.connect(self.get_temperature_Raman)
        # available_read_modes = ['FVB', 'Multi-track', 'Random track', 'Single track', 'Image']
        self.read_mode_Raman = 'Single track'
        self.acquisition_read_mode_Raman(self.read_mode_Raman)
        
        """ Prepare for file saving"""
        self.Datagroup_Raman = activeDatafile.create_group(name='AndorData', auto_increment=True)
        self.Datagroup_Raman.attrs.create('Version Date', str('Date of latest edit of Raman Rig2 = 15-09-2023'))
        
        """ update plot range"""
        self.update_figure_lim_Raman()
        self.SetLim_button_Raman.clicked.connect(self.update_figure_lim_Raman)
        

    """START: methods for use with OO spectrometer"""
    
    def set_default_int_times(self, default_t_in_ms):
        self.OOspectrometer.integration_time  = default_t_in_ms
        self.OOspectrometer.background_int_ref= default_t_in_ms
        self.OOspectrometer.background_int    = default_t_in_ms
        self.OOspectrometer.reference_int     = default_t_in_ms
        
    def save2txtfile(self, filename, integrationtime, wavelength, spectrum):
        now = datetime.now()
        current_date = now.strftime("%y%m%d")
        
        self.folder2save = os.path.join('C:\\Users\\Lab Di Martino\\Documents\\data\\' + self.OOspectrometer.CRSID)
        
        folders = []
        for root, dirs, files in os.walk(self.folder2save):
        
            for folder in dirs:
                if folder.startswith(current_date):
                    folders.append(folder)
        if not os.path.exists(self.folder2save + '\\' + str(folders[0])):
            os.makedirs(self.folder2save + '\\' + str(folders[0]))
        
        if spectrum.ndim == 1:
            # save single spectrum (no DF), supposed that they are for DF setting
            path2save = os.path.join(self.folder2save + '\\' + str(folders[0]) + '\\DF_bg_ref')
            fullfile = path2save + '\\DF_setting_' + filename
        elif spectrum.ndim > 1:
            path2save = os.path.join(self.folder2save + '\\' + str(folders[0]) + '\\DF_data')
            fullfile = path2save + '\\DF_' + filename
            
        if not os.path.exists(path2save):
            os.makedirs(path2save)
        
        # Save as increment file name
        i = 0        
        if os.path.exists(fullfile + '.txt'):
            while os.path.exists(fullfile + '_' + str(i) + '.txt'):
                i = i + 1
            print(i)
            outputfile = fullfile + '_' + str(i) + '.txt'
        else:
            outputfile = fullfile + '.txt'
            
        afile = open(outputfile, 'w')
        
        if spectrum.ndim == 1:
            # if single spectrum (no DF) will be saved
            saveData = np.vstack(np.transpose((\
                                               np.append([0], wavelength), \
                                               np.append([integrationtime], spectrum))))
        elif spectrum.ndim > 1:
            # if DF spectrum will be saved, all corresponding spectra will be
            # saved into the same file
            wavelength2save = np.append([0], wavelength)
            saveData = np.transpose(np.vstack((wavelength2save, spectrum)))
                
        np.savetxt(afile, saveData)
        afile.close()
    
    def RestoreDFsetting(self):
        now = datetime.now()
        current_date = now.strftime("%y%m%d")
        
        self.folder2save = os.path.join('C:\\Users\\Lab Di Martino\\Documents\\data\\' + self.OOspectrometer.CRSID)
        folders = []
        for root, dirs, files in os.walk(self.folder2save):
        
            for folder in dirs:
                if folder.startswith(current_date):
                    folders.append(folder)
                    
        path2save = os.path.join(self.folder2save + '\\' + str(folders[0]) + '\\DF_bg_ref')
        
        ''' Load background_ref'''
        # Always load the last saved spectra
        filename = 'background_ref'
        i = 0
        
        fullfile = path2save + '\\DF_setting_' + filename
        if os.path.exists(fullfile + '.txt'):
            outputfile = fullfile + '.txt'
            loadedfilei = ''
            while os.path.exists(fullfile + '_' + str(i) + '.txt'):
                i = i + 1
            if i > 0:
                i = i - 1
                outputfile = fullfile + '_' + str(i) + '.txt'
                loadedfilei = str(i)
            testI = np.loadtxt(outputfile)
                        
            load_x_data = []
            load_y_data = []
        
            for row in testI:
                load_x_data.append(row[0])
                load_y_data.append(row[1])
            
            self.OOspectrometer.background_int_ref = load_y_data[0].astype(float)
            self.OOspectrometer.background_ref = np.array(load_y_data[1:len(load_y_data)])
            print('Load background_ref ' + loadedfilei)
        else:
            print('No file found')
                
        ''' Load reference'''
        # Always load the last saved spectra
        filename = 'reference'
        i = 0
        
        fullfile = path2save + '\\DF_setting_' + filename
        if os.path.exists(fullfile + '.txt'):
            outputfile = fullfile + '.txt'
            loadedfilei = ''
            while os.path.exists(fullfile + '_' + str(i) + '.txt'):
                i = i + 1
            if i > 0:
                i = i - 1
                outputfile = fullfile + '_' + str(i) + '.txt'
                loadedfilei = str(i)
            testI = np.loadtxt(outputfile)
                        
            load_x_data = []
            load_y_data = []
        
            for row in testI:
                load_x_data.append(row[0])
                load_y_data.append(row[1])
            
            self.OOspectrometer.reference_int = load_y_data[0].astype(float)
            self.OOspectrometer.reference = np.array(load_y_data[1:len(load_y_data)])
            print('Load reference ' + loadedfilei)
        else:
            print('No file found')
        
        ''' Load background'''
        # Always load the last saved spectra
        filename = 'background'
        i = 0
        
        fullfile = path2save + '\\DF_setting_' + filename
        if os.path.exists(fullfile + '.txt'):
            outputfile = fullfile + '.txt'
            loadedfilei = ''
            while os.path.exists(fullfile + '_' + str(i) + '.txt'):
                i = i + 1
            if i > 0:
                i = i - 1
                outputfile = fullfile + '_' + str(i) + '.txt'
                loadedfilei = str(i)
            testI = np.loadtxt(outputfile)
                        
            load_x_data = []
            load_y_data = []
        
            for row in testI:
                load_x_data.append(row[0])
                load_y_data.append(row[1])
            
            self.OOspectrometer.background_int = load_y_data[0].astype(float)
            self.OOspectrometer.background = np.array(load_y_data[1:len(load_y_data)])
            print('Load background ' + loadedfilei)
        else:
            print('No file found')
        
    def read_bkg_for_ref_method(self):
        #use parallel thread to tell spectrometer to read spectrum. Mark spectrometer as in use until it's finished reading.
        take_spec_task = sub_take_spec(self.OOspectrometer)
        
        self.OOspectrometer.OO_occupied = 1
        self.pool.start(take_spec_task)
        # Wait till spectrum is taken
        while self.OOspectrometer.OO_occupied > 0 :
            pass
        #once finished:
        if self.OOspectrometer.OO_occupied == 0:
            #store the spectrum
            self.OOspectrometer.background_ref = self.remove_invalid_values_from_spectrum(take_spec_task.current_spec)
            self.OOspectrometer.background_int_ref = self.OOspectrometer.integration_time
            self.background_subtracted_ref.setCheckState(QtCore.Qt.Checked)
            ''' Trung Feb 2024'''
            # Save into text file
            self.save2txtfile('background_ref', \
                              self.OOspectrometer.background_int_ref, \
                              self.wl_in_nm, \
                              self.OOspectrometer.background_ref)
            print('background for reference has been stored')
            
    def clear_bkg_for_ref_method(self):
        self.OOspectrometer.background_ref = np.zeros(len(self.wl_in_nm ))
        self.background_subtracted_ref.setCheckState(QtCore.Qt.Unchecked)
    
    def read_ref_method(self):
        #use parallel thread to tell spectrometer to read spectrum. Mark spectrometer as in use until it's finished reading.
        take_spec_task = sub_take_spec(self.OOspectrometer)
       
        self.OOspectrometer.OO_occupied = 1
        self.pool.start(take_spec_task)
        # Wait till spectrum is taken
        while self.OOspectrometer.OO_occupied > 0:
            pass
        #once finished:
        if self.OOspectrometer.OO_occupied == 0:
            self.OOspectrometer.reference  = self.remove_invalid_values_from_spectrum(take_spec_task.current_spec)
            self.OOspectrometer.reference_int = self.OOspectrometer.integration_time
            self.referenced.setCheckState(QtCore.Qt.Checked)
            ''' Trung Feb 2024'''
            # Save into text file
            self.save2txtfile('reference', \
                              self.OOspectrometer.reference_int, \
                              self.wl_in_nm, \
                              self.OOspectrometer.reference)
            print('reference has been stored')


    def clear_ref_method(self):
        self.OOspectrometer.reference = np.ones(len(self.wl_in_nm ))
        self.referenced.setCheckState(QtCore.Qt.Unchecked)
        
    def read_background_method(self):
        #use parallel thread to tell spectrometer to read spectrum. Mark spectrometer as in use until it's finished reading.
        take_spec_task = sub_take_spec(self.OOspectrometer)
        
        self.OOspectrometer.OO_occupied = 1
        self.pool.start(take_spec_task)
        # Wait till spectrum is taken
        while self.OOspectrometer.OO_occupied > 0:
            pass
        #once finished:
        if self.OOspectrometer.OO_occupied == 0:
            #store the spectrum
            self.OOspectrometer.background  = self.remove_invalid_values_from_spectrum(take_spec_task.current_spec)
            self.OOspectrometer.background_int= self.OOspectrometer.integration_time
            self.background_subtracted.setCheckState(QtCore.Qt.Checked)
            ''' Trung Feb 2024'''
            # Save into text file
            self.save2txtfile('background', \
                              self.OOspectrometer.background_int, \
                              self.wl_in_nm, \
                              self.OOspectrometer.background)
            print('background has been stored')

        
    def clear_background_method(self):
        self.OOspectrometer.background = np.zeros(len(self.wl_in_nm ))
        self.background_subtracted.setCheckState(QtCore.Qt.Unchecked)

    def take_spec(self):    

        # read the current field strength and sample temperature
        if self.MiPS_handle != None:
            current_field = sub_read_current_field(self.MiPS_handle, self.cFieldT)
            self.pool.start(current_field)
        if self.MiTC_handle  != None:
            current_temp = sub_read_current_sample_temp(self.MiTC_handle, self.cSampleTemp)
            self.pool.start(current_temp)
        
        #hold cryostat stage steady during integration of signal by the spectrometer
        if self.stage_enabled == 'Cryostat stage':
            # change the movemode of channel 1 to absolute and read the current stage position 
            move_mode = smaract_package.MoveMode.CL_ABSOLUTE
            self.smaract_set_movemode(1, move_mode)
            self.current_z_position = smaract_package.GetProperty_i64(self.Smaract_stage, 1, smaract_package.Property.POSITION)
            stage_on_hold = sub_stage_on_hold(self.OOspectrometer, self.Smaract_stage, self.current_z_position, 1)
            self.pool.start(stage_on_hold)
        #if using olympus stage just store curret z position of stage.
        if self.stage_enabled == 'Olympus stage':
            global stage_in_use
            while stage_in_use ==True:            
                pass
            if stage_in_use ==False:    
                self.current_z_position = self.SMC_stage.get_position(3)
        
        # Use helper class to setup task to take the spectrum in the background
        take_spec_task = sub_take_spec(self.OOspectrometer)
        # mark spectrometer in use, and set task running
        
        self.OOspectrometer.OO_occupied = 1
        self.pool.start(take_spec_task)
        
        # Wait till spectrum is taken
        while self.OOspectrometer.OO_occupied > 0:
            pass
        #once finished:
        if self.OOspectrometer.OO_occupied == 0:
            #store the spectrum
            self.current_raw_spec = None
            self.current_processed_spec = None
            self.current_raw_spec = self.remove_invalid_values_from_spectrum(take_spec_task.current_spec)
            self.current_processed_spec = self.process_spec(self.current_raw_spec)
            #change smaract z-axis movemode back to the normal 'step' setting
            if self.stage_enabled == 'Cryostat stage':
                channel = 1
                self.smaract_set_movemode(channel, self.smaract_movemode)
            #store solenoid params, if applicable
            if self.MiPS_handle != None:
                self.current_field = current_field.cFieldT_value
            elif self.MiPS_handle == None:
                self.current_field =0
            if self.MiTC_handle != None:  
                self.current_temp = current_temp.cSampleTemp_value
            elif self.MiTC_handle == None: 
                self.current_temp  = 0

    def remove_invalid_values_from_spectrum(self,spectrum):
        inf_truth = np.isinf(spectrum)
        nan_truth = np.isnan(spectrum)
        i=0
        while i<len(spectrum):
            if inf_truth[i]==True or nan_truth[i]==True:
                spectrum[i]=0
            i=i+1
        return spectrum
            
    def take_spec_and_update_display(self):
        # take_spec method allows parallel processing
        self.take_spec()
        #this line in twice to get rid of weird timing glitch...temporary patch...
        #self.take_spec()
        
        while self.OOspectrometer.OO_occupied > 0:
            pass
        if self.OOspectrometer.OO_occupied == 0:
            self.Spectra_display.canvas.ax.cla()
            p = self.Spectra_display.canvas.ax.plot(self.wl_in_nm, self.current_processed_spec)
            self.Spectra_display.canvas.ax.set_ylabel('DF signal', fontsize = 16)
            self.Spectra_display.canvas.ax.set_xlabel('Wavelength (nm)', fontsize = 16)
            self.Spectra_display.canvas.ax.tick_params(axis = "y", direction = "in")
            self.Spectra_display.canvas.ax.tick_params(axis = "x", direction = "in", top = False)
            
            minY = min(float(self.YLim1_box.value()), float(self.YLim2_box.value()))
            maxY = max(float(self.YLim1_box.value()), float(self.YLim2_box.value()))
            if minY == maxY:
                pos1 = min(np.argmin(self.wl_in_nm - 500), np.argmin(self.wl_in_nm - 800))
                pos2 = max(np.argmin(self.wl_in_nm - 500), np.argmin(self.wl_in_nm - 800))
                minY = min(self.current_processed_spec[pos1:pos2 + 1])
                maxY = max(self.current_processed_spec[pos1:pos2 + 1])
            
            self.Spectra_display.canvas.ax.set_ylim([minY, maxY])
            minX = min(float(self.XLim1_box.value()), float(self.XLim2_box.value()))
            maxX = max(float(self.XLim1_box.value()), float(self.XLim2_box.value()))
            self.Spectra_display.canvas.ax.set_xlim([minX, maxX])
            self.Spectra_display.canvas.fig.tight_layout()
            self.Spectra_display.canvas.draw()
#            self.spec_plot.clear()
#            self.spec_plot.plot(self.wl_in_nm, self.current_processed_spec)
            print('single spec display has been updated')
        
    def SetfigRange(self):
        minY = min(float(self.YLim1_box.value()), float(self.YLim2_box.value()))
        maxY = max(float(self.YLim1_box.value()), float(self.YLim2_box.value()))
        if minY == maxY:
            pos1 = min(np.argmin(self.wl_in_nm - 500), np.argmin(self.wl_in_nm - 800))
            pos2 = max(np.argmin(self.wl_in_nm - 500), np.argmin(self.wl_in_nm - 800))
            minY = min(self.current_processed_spec[pos1:pos2 + 1])
            maxY = max(self.current_processed_spec[pos1:pos2 + 1])
        
        self.Spectra_display.canvas.ax.set_ylim([minY, maxY])
        minX = min(float(self.XLim1_box.value()), float(self.XLim2_box.value()))
        maxX = max(float(self.XLim1_box.value()), float(self.XLim2_box.value()))
        self.Spectra_display.canvas.ax.set_xlim([minX, maxX])
        self.Spectra_display.canvas.fig.tight_layout()
        self.Spectra_display.canvas.draw()
        
        self.SpectraStack_display.canvas.ax.set_ylim([minY, maxY])
        minX = min(float(self.XLim1_box.value()), float(self.XLim2_box.value()))
        maxX = max(float(self.XLim1_box.value()), float(self.XLim2_box.value()))
        self.SpectraStack_display.canvas.ax.set_xlim([minX, maxX])
        self.SpectraStack_display.canvas.fig.tight_layout()
        self.SpectraStack_display.canvas.draw()
        
    def process_spec(self, a_spectrum):

        #scale bkg/ref spectra according to ratio of integration times with data spectrum being processed
        bkg_for_ref_scaling_factor = self.OOspectrometer.integration_time/self.OOspectrometer.background_int_ref
        bkg_scaling_factor  =  (self.OOspectrometer.integration_time/self.OOspectrometer.background_int)
        ref_scaling_factor  = (self.OOspectrometer.integration_time/self.OOspectrometer.reference_int)
        
        #if no bkg for ref taken, just use zeros
        if self.background_subtracted_ref.isChecked()==False:
            print('note - no bkg for ref has been stored')
            self.OOspectrometer.background_ref = np.zeros(len(a_spectrum))
        #if no reference taken, just used ones
        if self.referenced.isChecked()==False:
            print('note - no ref has been stored')
            ref_with_bkg_subtracted = np.ones(len(a_spectrum))
        #if reference has been taken, subtract bkg for ref (both quantities scaled according to integration times)    
        elif self.referenced.isChecked()==True:
            ref_with_bkg_subtracted  = (self.OOspectrometer.reference * ref_scaling_factor) - (self.OOspectrometer.background_ref * bkg_for_ref_scaling_factor)
            #check whether white or metal reference is used - use/don't use ratio curve accordingly
            if self.referencing_type == 'Gold':
                ref_with_bkg_subtracted  = ref_with_bkg_subtracted * self.ratio_curve_values 

        #calculate processed spectrum
        if self.background_subtracted.isChecked()==False:
            print('note - no bkg has been stored')
            scaled_bkg = np.zeros(len(a_spectrum))
        elif self.background_subtracted.isChecked()==True:
            scaled_bkg = self.OOspectrometer.background * bkg_scaling_factor
        processed_spectrum = (a_spectrum - scaled_bkg)/(ref_with_bkg_subtracted)

        return processed_spectrum 

    """connect spectrometer control and readout buttons"""
#    def gui_set_tec_temperature(self):
#        self.OOspectrometer.tec_temperature = float(self.set_tec_temperature_LineEdit.text().strip())
        
    def gui_read_tec_tempeature(self):
        self.tec_temperature_lcdNumber.display(float(self.OOspectrometer.tec_temperature)) 
        
    def gui_set_integration_time(self):
        self.OOspectrometer.integration_time = float(self.set_integration_time_LineEdit.text().strip())
        self.OOspectrometer.set_integration_time(self.OOspectrometer.integration_time)
        #introduce wait with aim of avoiding weird glitch with spectrum update
        time.sleep(self.OOspectrometer.integration_time/1000)
        print(' ')
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
        
        """ Save into txt file """
        filename = str(self.description_LineEdit.text().strip())
        if filename == 'Description...':
            filename == 'NP' + str(self.current_temp) + 'K_'+ \
                    str(self.current_field) + 'T_' + \
                    str(self.current_z_position) + 'mm'
        
        # spectra of processed, NP, background, reference, bgref
        process2save = np.append([self.OOspectrometer.integration_time], self.current_processed_spec)
        NP2save = np.append([self.OOspectrometer.integration_time], self.current_raw_spec)
        background2save = np.append([self.OOspectrometer.background_int], self.OOspectrometer.background)
        reference2save = np.append([self.OOspectrometer.reference_int], self.OOspectrometer.reference)
        BGreference2save = np.append([self.OOspectrometer.background_int_ref], self.OOspectrometer.background_ref)
        
        spec2save = np.array(np.vstack((process2save, NP2save, background2save, reference2save, BGreference2save)))

        self.save2txtfile(filename, \
                              0, \
                              self.wl_in_nm, \
                              spec2save)
        
        print(' ')
        print('A single DF spectrum has been saved.' )
        print('Solenoid params =  '+str(self.current_field) + ' T and ' + str(self.current_temp) + ' K.')        
        
        
        
        
        
    def t_series_method(self):
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
        print(' ')
        print('  BEGIN TIME SERIES')
        print('t_series_delay = '+str(t_series_delay*1000)+' ms.')
        print(' self.OOspectrometer.integration_time = '+str(self.OOspectrometer.integration_time)+' ms.')
        print(' - - -')
        while spectrum_counter < t_series_num_spec:
            print('spectrum_counter   = '+str(spectrum_counter ))
            print('time_counter   = '+str(time_counter ))
            self.take_spec()
            self.t_series_times             = self.t_series_times             + [time_counter ]
            self.raw_t_series_spectra       = self.raw_t_series_spectra       + [self.current_raw_spec]
            self.processed_t_series_spectra = self.processed_t_series_spectra + [self.current_processed_spec]
            time.sleep(t_series_delay)
            time_counter = time_counter + self.OOspectrometer.integration_time + (t_series_delay*1000)
            spectrum_counter  = spectrum_counter + 1
            print(' - ')

        activeTSeriesDataset = self.time_series_group.create_dataset('SpectrumTimeSeries_%d', data = self.raw_t_series_spectra)
        
        activeTSeriesDataset.attrs.create('t-series num spectra', t_series_num_spec)
        activeTSeriesDataset.attrs.create('t-series spec end-to-start delay', t_series_delay)
        activeTSeriesDataset.attrs.create("Wavelengths", self.wl_in_nm)
        activeTSeriesDataset.attrs.create('Integration_Time', self.OOspectrometer.integration_time)         
        activeTSeriesDataset.attrs.create("z-position", self.current_z_position)
        #activeTSeriesDataset.attrs.create("Processed Spectra", self.processed_t_series_spectra )
        
        self.read_spec_description()
        activeTSeriesDataset.attrs.create("Description", self.current_spec_description)
        activeTSeriesDataset.attrs.create('bkg_for_ref_int_t', self.OOspectrometer.background_int_ref)
        activeTSeriesDataset.attrs.create('bkg_int_t', self.OOspectrometer.background_int)
        activeTSeriesDataset.attrs.create('ref_int_t', self.OOspectrometer.reference_int)
        activeTSeriesDataset.attrs.create('bkg_for_ref', self.OOspectrometer.background_ref )
        activeTSeriesDataset.attrs.create('background', self.OOspectrometer.background)
        activeTSeriesDataset.attrs.create('reference', self.OOspectrometer.reference)
        
        activeTSeriesDataset.attrs.create('solenoid temperature', self.current_temp)
        activeTSeriesDataset.attrs.create('solenoid field', self.current_field)
        print('  END TIME SERIES')
        print(' ')
        
    """END: methods for use with OO spectrometer"""
            
    """START: methods for z-stacks using OO spectrometer"""
    def z_stack_method(self):
        self.Spectra_display.canvas.ax.cla()
#        self.stack_plot.clear()
        
        self.stack_step_size = int(self.zstep_box_2.value())
        self.stack_num_steps = int(self.num_steps_box.value())
        self.stack_z_values= []
        
        global z_stack_running
        z_stack_running = True
                
        self.processed_stack_spectra = []
        self.raw_stack_spectra = []
        step_num = 0
        print(' ' )
        print('BEGIN Z STACK')
        while step_num < self.stack_num_steps:
            print('step_num = '+str(step_num))
            #take_spec() method uses sub_take_spec class to run parallel thread for read_spectrum() instruction to the spectrometer
            self.take_spec()
            #after taking 1st spec: store stage z as initial position
            if step_num==0:
                initial_z = self.current_z_position
            #for all steps, store z_position 
            self.stack_z_values = self.stack_z_values + [self.current_z_position]
            #check OO_occupied--> will change from True to False once spectrometer is finished reading the spectrum

            while self.OOspectrometer.OO_occupied > 0:
                pass
            if self.OOspectrometer.OO_occupied == 0:
                #once spectrum finished, move to new z for next spectrum (unless on final spec)
                if step_num==self.stack_num_steps:
                    pass
                elif step_num!=self.stack_num_steps:
                    if self.stage_enabled == 'Cryostat stage':
                        self.smaract_move(1, self.stack_step_size*(-1))
                    elif self.stage_enabled == 'Olympus stage':
                        self.SMC100_move_zaxis(self.stack_step_size)
                #save information gathered on this z_increment
                self.processed_stack_spectra = self.processed_stack_spectra + [self.current_processed_spec]
                self.raw_stack_spectra = self.raw_stack_spectra + [self.current_raw_spec]
                p = self.SpectraStack_display.canvas.ax.plot(self.wl_in_nm, self.current_processed_spec)
                self.SpectraStack_display.canvas.ax.set_ylabel('DF signal', fontsize = 16)
                self.SpectraStack_display.canvas.ax.set_xlabel('Wavelength (nm)', fontsize = 16)
                self.SpectraStack_display.canvas.ax.tick_params(axis = "y", direction = "in")
                self.SpectraStack_display.canvas.ax.tick_params(axis = "x", direction = "in", top = False)
                
                minY = min(float(self.YLim1_box.value()), float(self.YLim2_box.value()))
                maxY = max(float(self.YLim1_box.value()), float(self.YLim2_box.value()))
                if minY == maxY:
                    pos1 = min(np.argmin(self.wl_in_nm - 500), np.argmin(self.wl_in_nm - 800))
                    pos2 = max(np.argmin(self.wl_in_nm - 500), np.argmin(self.wl_in_nm - 800))
                    minY = min(self.current_processed_spec[pos1:pos2 + 1])
                    maxY = max(self.current_processed_spec[pos1:pos2 + 1])
                    
                self.SpectraStack_display.canvas.ax.set_ylim([minY, maxY])
                minX = min(float(self.XLim1_box.value()), float(self.XLim2_box.value()))
                maxX = max(float(self.XLim1_box.value()), float(self.XLim2_box.value()))
                self.SpectraStack_display.canvas.ax.set_xlim([minX, maxX])
                self.Spectra_display.canvas.fig.tight_layout()
                self.SpectraStack_display.canvas.draw()
#                self.stack_plot.plot(self.wl_in_nm, self.current_processed_spec)
                
                #increment step counter
                step_num = step_num+1
        print('END Z STACK')
        print(' ' )
        final_z = self.current_z_position
        
        #save the raw spectra as dataset, and save other relevant info as attributes to the dataset
        activeStackDataset = self.Stacked_Spectra_group.create_dataset('stackedDFSpectra_%d', data = self.raw_stack_spectra)
        activeStackDataset.attrs.create('stack step size', self.stack_step_size)
        activeStackDataset.attrs.create('stack step number', self.stack_num_steps)
        activeStackDataset.attrs.create('stack z values', self.stack_z_values )
        activeStackDataset.attrs.create("Initial z position", initial_z)
        activeStackDataset.attrs.create("Final z position", final_z)
        
        self.read_spec_description()
        activeStackDataset.attrs.create("Description", self.current_spec_description)
        activeStackDataset.attrs.create("Wavelengths", self.wl_in_nm)
        activeStackDataset.attrs.create('Integration_Time', self.OOspectrometer.integration_time)
        #activeStackDataset.attrs.create("Processed Spectra", self.processed_stack_spectra )
        
        activeStackDataset.attrs.create('bkg_for_ref_int_t', self.OOspectrometer.background_int_ref)
        activeStackDataset.attrs.create('bkg_int_t', self.OOspectrometer.background_int)
        activeStackDataset.attrs.create('ref_int_t', self.OOspectrometer.reference_int)
        activeStackDataset.attrs.create('bkg_for_ref', self.OOspectrometer.background_ref )
        activeStackDataset.attrs.create('background', self.OOspectrometer.background)
        activeStackDataset.attrs.create('reference', self.OOspectrometer.reference)
        
        activeStackDataset.attrs.create('solenoid temperature', self.current_temp)
        activeStackDataset.attrs.create('solenoid field', self.current_field)
        
        self.make_stack_plot()

            
    def make_stack_plot(self):
        spec_nums = np.arange(0, self.stack_num_steps, 1)
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
            print('Olympus stage disconnected.')
            
    def Set_New_Origin(self):
        #this will take current coordinates of the three smaracts and save them as 'origin' - note that correct frame of reference must be used
        if self.stage_enabled == 'Cryostat stage':
            print(' Saving coords as origin - check frame of reference is use is the desired one.')
            z_position = smaract_package.GetProperty_i64(self.Smaract_stage, 1, smaract_package.Property.POSITION)
            left_right_position = smaract_package.GetProperty_i64(self.Smaract_stage, self.smaract_channel_left, smaract_package.Property.POSITION)
            up_down_position = smaract_package.GetProperty_i64(self.Smaract_stage, self.smaract_channel_up, smaract_package.Property.POSITION)
            
            self.chosen_coord_origin={}
            self.chosen_coord_origin['z_position']          = z_position
            self.chosen_coord_origin['left_right_position'] = left_right_position
            self.chosen_coord_origin['up_down_position']    = up_down_position
            
            print('origin coords stored as:')
            print('left/right axis = '+ str(self.chosen_coord_origin['z_position']/10**(6))         +' um.' )
            print('up/down axis    = '+ str(self.chosen_coord_origin['left_right_position']/10**(6))+' um.' )
            print('z-axis          = '+ str(self.chosen_coord_origin['up_down_position']/10**(6))   +' um.' )
            print('   -   -   -  ')
            
            self.origin_stored_checkbox.setChecked(True)
        
        elif self.stage_enabled == 'Olympus stage':
            print('origin not yet implemented for SMC100')
            z_position = self.SMC_stage.get_position(3)
            z_position = z_position[0]

        
    def Go_To_Stage_Origin(self):
        if self.stage_enabled == 'Cryostat stage':
            
            print('moving to origin cryostat stage - backlash calibration needed?')
            
            z_position = smaract_package.GetProperty_i64(self.Smaract_stage, 1, smaract_package.Property.POSITION)
            left_right_position = smaract_package.GetProperty_i64(self.Smaract_stage, self.smaract_channel_left, smaract_package.Property.POSITION)
            up_down_position = smaract_package.GetProperty_i64(self.Smaract_stage, self.smaract_channel_up, smaract_package.Property.POSITION)

            required_z_steps   = int( (self.chosen_coord_origin['z_position']/10**(6))          - ( z_position/10**(6) ) )
            required_l_r_steps = int( (self.chosen_coord_origin['left_right_position']/10**(6)) - ( left_right_position /10**(6) ) )
            required_u_d_steps = int( (self.chosen_coord_origin['up_down_position']/10**(6))    - ( up_down_position/10**(6) ) )
 
            self.smaract_move(channel=self.smaract_channel_left, move_value = required_l_r_steps)
            self.smaract_move(channel=1,                         move_value = required_z_steps)            
            self.smaract_move(channel=self.smaract_channel_up,   move_value = required_u_d_steps)
            print(' ')
            
        elif self.stage_enabled == 'Olympus stage':
            print('origin not yet implemented for SMC100')
            self.SMC_stage.move(self.current_z_setpoint, '3', relative=False)
            print('SMC stage returned to z-axis setpoint')
            print('Position axis 3: ' + str(self.current_z_setpoint) + ' mm')
            

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
        self.print_stage_positions()
        
    def smaract_moveright(self):
        # channel 0 (-)
        self.smaract_move_value = int(self.xystep_Box.value())
        self.smaract_move(self.smaract_channel_right, self.smaract_move_value*self.smaract_direction_right)
        self.print_stage_positions()
        
    def smaract_moveup(self):
        # channel 2 (-)
        self.smaract_move_value = int(self.xystep_Box.value())
        self.smaract_move(self.smaract_channel_up, self.smaract_move_value*self.smaract_direction_up)
        self.print_stage_positions()
        
    def smaract_movedown(self):
        # channel 2 (+)
        self.smaract_move_value = int(self.xystep_Box.value())
        self.smaract_move(self.smaract_channel_down, self.smaract_move_value*self.smaract_direction_down)
        self.print_stage_positions()
        
    def smaract_movetoward(self):
        # channel 1 (+)
        self.smaract_move_value = int(self.zstep_box.value())
        self.smaract_move(1, self.smaract_move_value)
        self.print_stage_positions()
        
    def smaract_moveaway(self):
        # channel 1 (-)
        self.smaract_move_value = int(self.zstep_box.value())
        self.smaract_move(1, self.smaract_move_value*(-1))
        self.print_stage_positions()
        
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
        global stage_in_use
        stage_in_use =True
        self.pool.start(do_initialise) #once initialisation is finished, stage_in_use is set to False
        self.stage_feature_button.setEnabled(False)
        
    def SMC100_move_left(self):
        axis = 1
        goto_value = float(self.xystep_Box.value())*100/1e6
        global stage_in_use
        while stage_in_use==True:
            pass
        if stage_in_use==False:
            do_move = sub_SMC100_move(self.SMC_stage, str(axis), goto_value)
            stage_in_use =True
            self.pool.start(do_move)
    
    def SMC100_move_right(self):
        axis = 1
        goto_value = -float(self.xystep_Box.value())*100/1e6
        global stage_in_use
        while stage_in_use==True:
            pass
        if stage_in_use==False:
            do_move = sub_SMC100_move(self.SMC_stage, str(axis), goto_value)
            stage_in_use =True
            self.pool.start(do_move)
    
    def SMC100_move_up(self):
        axis = 2
        goto_value = float(self.xystep_Box.value())*100/1e6
        global stage_in_use
        while stage_in_use==True:
            pass
        if stage_in_use==False:
            do_move = sub_SMC100_move(self.SMC_stage, str(axis), goto_value)
            stage_in_use =True
            self.pool.start(do_move)
    
    def SMC100_move_down(self):
        axis = 2
        goto_value = -float(self.xystep_Box.value())*100/1e6
        global stage_in_use
        while stage_in_use==True:
            pass
        if stage_in_use==False:
            do_move = sub_SMC100_move(self.SMC_stage, str(axis), goto_value)
            stage_in_use =True
            self.pool.start(do_move)

    def SMC100_move_towardslens(self):
        axis = 3
        goto_value = float(self.zstep_box.value())*50/1e6
        global stage_in_use
        while stage_in_use==True:
            pass
        if stage_in_use==False:
            do_move = sub_SMC100_move(self.SMC_stage, str(axis), goto_value)
            stage_in_use =True
            self.pool.start(do_move)
    
    def SMC100_move_awaylens(self):
        axis = 3
        goto_value = -float(self.zstep_box.value())*50/1e6
        global stage_in_use
        while stage_in_use==True:
            pass
        if stage_in_use==False:
            do_move = sub_SMC100_move(self.SMC_stage, str(axis), goto_value)
            stage_in_use =True
            self.pool.start(do_move)
        
    def SMC100_move_zaxis(self, step_size):
        axis = 3
        goto_value = -float(step_size)*50/1e6
        global stage_in_use
        while stage_in_use==True:
            pass
        if stage_in_use==False:
            do_move = sub_SMC100_move(self.SMC_stage, str(axis), goto_value)
            stage_in_use =True
            self.pool.start(do_move)
        
    def SMC100_move_mid(self):
        if self.stage_enabled == 'Cryostat stage':
            print("no 'go to midpoint' function defined for cryostat stage")
        elif self.stage_enabled == 'Olympus stage':
            do_move_mid = sub_SMC100_move_mid(self.SMC_stage)
            global stage_in_use
            stage_in_use =True
            self.pool.start(do_move_mid)
        
    def print_stage_positions(self):
        #note - cryostat stage case is untested
        #note - need to add units to position values
        
        if self.stage_enabled == 'Cryostat stage':
            z_position = smaract_package.GetProperty_i64(self.Smaract_stage, 1, smaract_package.Property.POSITION)
            left_right_position = smaract_package.GetProperty_i64(self.Smaract_stage, self.smaract_channel_left, smaract_package.Property.POSITION)
            up_down_position = smaract_package.GetProperty_i64(self.Smaract_stage, self.smaract_channel_up, smaract_package.Property.POSITION)
            self.z_pos_display.display(float(z_position/10**(6)))
            self.lr_pos_display.display(float(left_right_position/10**(6)))
            self.ud_pos_display.display(float(up_down_position/10**(6)))
            # 500000 units on smaract axis = 500nm
            # 1000 units on smaract axis = 1nm
            # axis/1000 = nm = 10^-3 um
            # 1nm
            
            button_text = self.stage_feature_button.text()
            
            
            print('     - - - - -     ')
            print(button_text)
            print('left/right axis = '+str(left_right_position/10**(6))+' um (moves in ?nm increments)')
            print('up/down axis    = '+str(up_down_position/10**(6))+' um (moves in ?nm increments)')
            print('z-axis          ='+str(z_position/10**(6))+' um (moves in ?nm increments)')
            print('     - - - - -     ')
            
        elif self.stage_enabled == 'Olympus stage':
            global stage_in_use
            while stage_in_use==True:
                pass
            if stage_in_use==False:
                z_position = self.SMC_stage.get_position(3)
                left_right_position = self.SMC_stage.get_position(1)
                up_down_position = self.SMC_stage.get_position(2)
                print('          ')
                print('Olympus stage positions:')
                print('     - - - - -     ')
                print('left/right axis = '+str(left_right_position)+'mm (moves in 100nm increments)')
                print('up/down axis = '+str(up_down_position)+'mm (moves in 100nm increments)')
                print('z-axis = '+str(z_position)+'mm (moves in 50nm increments)')
                print('     - - - - -     ')
                self.z_pos_display.display(float(z_position[0]))
                self.lr_pos_display.display(float(left_right_position[0]))
                self.ud_pos_display.display(float(up_down_position[0]))
        
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

    """ Raman part, Trung 11.04.2024 ------------------------------------------ """
    
    def Raman_scan_direction(self):
        available_modes = ['Direction Left', 'Direction Up', 'Direction Right', 'Direction Down']
        button_text = self.Switch_scan_direction_button_Raman.text()
        if available_modes.index(button_text) == len(available_modes)-1:
            nextpos = 0
        else:
            nextpos = available_modes.index(button_text) + 1
        self.Switch_scan_direction_button_Raman.setText(available_modes[nextpos])
    
    def Capture_xy_scan_Raman(self):
        available_modes = ['Direction Left', 'Direction Up', 'Direction Right', 'Direction Down']
        button_text = self.Switch_scan_direction_button_Raman.text()
        nrstep = self.xy_scan_step_button_Raman.value()
        
        spectra_xy_scan = []
        
        """ Capture spectrum """
        
        if self.Andor.Andor_occupied > 0:
            return
        
        self.Camerastate_output_Raman.setText('Capturing ...')
        self.Camerastate_output_Raman.setStyleSheet('background-color: rgb(255, 140, 140);')
        
        currentstep = 0
        while currentstep <= nrstep:
            # Take spectrum
            imageArray, num_of_images, image_shape = self.Andor.capture()
            imageArray.reverse()
            
            spectrum_Raman = np.array(imageArray)
            # add into array
            spectra_xy_scan.append(spectrum_Raman.tolist())
            
            # move stage to the next position
            
            if button_text == available_modes[0]:
                self.SMC100_move_left()
            elif button_text == available_modes[1]:
                self.SMC100_move_up()
            elif button_text == available_modes[2]:
                self.SMC100_move_right()
            elif button_text == available_modes[3]:
                self.SMC100_move_down()
            while stage_in_use == True:
                pass
            currentstep = currentstep + 1

        # plot the spectra 
        self.Spectra_display_Raman.canvas.ax.cla()

        for i in range(0, len(spectra_xy_scan)):
            self.Spectra_display_Raman.canvas.ax.plot(wavelength_Raman, np.array(spectra_xy_scan[i]), label = str(i))
        self.Spectra_display_Raman.canvas.ax.legend()
        self.update_figure_lim_Raman()
        
        self.Camerastate_output_Raman.setText('Done!')
        self.Camerastate_output_Raman.setStyleSheet('background-color: rgb(180, 227, 255);')
    
    # Spectrometer parts
    def get_x_axis(self):
        return self.kymera.GetCalibration()[::-1]
    
    def set_centralwavelength_Raman(self):
        self.kymera.SetWavelength(self.Centralwavelength_box_Raman.value())
        global wavelength_Raman
        wavelength_Raman = np.array(self.get_x_axis())
        self.kymera.wavelength = wavelength_Raman
        self.update_figure_lim_Raman()
        
    # Camera parts

    def Andor_occupied_reset(self):
        global Abort_live_mode
        Abort_live_mode = 0
        
    def get_temperature_Raman(self):
       
        if self.Andor.Andor_occupied > 0:
            return
        
        read_temp = sub_get_temperature_Raman(self.Andor_temp_button_Raman, self.Andor)
        
        self.Andor.Andor_occupied = 0
        self.pool.start(read_temp)

    def acquisition_read_mode_Raman(self, readmode):
        if self.Andor.Andor_occupied > 0:
            print('No change applied')
            return
        
        available_modes = ['Single', 'Accumulate', 'Kinetic', 'Fast Kinetic']
        currentMode = 'Single'
        self.Andor.set_andor_parameter('AcquisitionMode', available_modes.index(currentMode) + 1)
        print('Acquisition mode set to Single')
        
        available_modes = ['FVB', 'Multi-track', 'Random track', 'Single track', 'Image']
        currentMode = readmode
        self.Andor.set_andor_parameter('ReadMode', available_modes.index(currentMode))
        print('Read mode set to ' + currentMode)
        
    def Switch_read_mode_Raman(self):
        
        if self.Andor.Andor_occupied > 0:
            print('No change applied')
            return
        
        button_text = self.Readmode_button_Raman.text()
        if button_text == 'Image':
            self.read_mode_Raman = 'Single track'
            self.Readmode_button_Raman.setText('Single')
            self.acquisition_read_mode_Raman(self.read_mode_Raman)
        elif button_text == 'Single':
            self.read_mode_Raman = 'Image'
            self.Readmode_button_Raman.setText('Image')
            self.acquisition_read_mode_Raman(self.read_mode_Raman)
            
    def set_exposureTime_Raman(self):
        if self.Andor.Andor_occupied > 0:
            print('No change applied')
            return
        
        self.Andor.set_andor_parameter('Exposure', float(self.Exposure_box_Raman.value()))
        print('Exposure time set to ' + str(round(self.Andor.Exposure, 2)) + ' s')
        self.Exposure_box_Raman.setValue(round(self.Andor.Exposure, 2))
        
    def CaptureSpectrum_Raman(self):

        if self.Andor.Andor_occupied > 0:
            return

        self.Camerastate_output_Raman.setText('Capturing ...')
        
        if self.read_mode_Raman == 'Single track':
            start_capture = sub_capture_Raman(self.Andor, self.Spectra_display_Raman, \
                                    self.Camerastate_output_Raman, \
                                    self.YLim1_box_Raman, self.YLim2_box_Raman)
        elif self.read_mode_Raman == 'Image':
            start_capture = sub_captureImage_Raman(self.Andor, self.Spectra_display_Raman, \
                                    self.Camerastate_output_Raman, \
                                    self.YLim1_box_Raman, self.YLim2_box_Raman)
        self.Andor.Andor_occupied = 1
        self.pool.start(start_capture)
        
    def Capture_and_save_Raman(self):
        if self.Andor.Andor_occupied > 0:
            return
        
        self.Camerastate_output_Raman.setText('Capturing ...')
        """ save data into text file"""
        now = datetime.now()
        current_date = now.strftime("%y%m%d")
        
        self.folder2save = os.path.join('C:\\Users\\Lab Di Martino\\Documents\\data\\' + self.Andor.CRSID)
        
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
            
        
        start_capture = sub_capture_and_save_Raman(self.Datagroup_Raman, self.Andor, self.Spectra_display_Raman, \
                                             self.Camerastate_output_Raman, \
                                             self.Repeat_box_Raman.value(), self.Filename_text_Raman, \
                                             self.Description_text_Raman, self.YLim1_box_Raman, \
                                             self.YLim2_box_Raman, path2save)
        self.Andor.Andor_occupied = 1
        self.pool.start(start_capture)
        
    def Abort_Raman(self):
       
        global Abort_live_mode
        Abort_live_mode = 1
        
    def Live_mode_Raman(self):
        
        if self.Andor.Andor_occupied > 0:
            return
        
        self.Camerastate_output_Raman.setText('Capturing ...')
        self.Camerastate_output_Raman.setStyleSheet('background-color: rgb(255, 140, 140);')
        start_capture = sub_Live_mode_Raman(self.Andor, self.Spectra_display_Raman, \
                                             self.Camerastate_output_Raman,
                                             self.YLim1_box_Raman, self.YLim2_box_Raman)
        self.Andor.Andor_occupied = 1
        self.pool.start(start_capture)
        
    
    def setROI_Raman(self):
        
        if self.Andor.Andor_occupied > 0:
            print('No change applied')
            return
        
        central_row = self.Centralrow_box_Raman.value()
        rows = self.Rows_box_Raman.value()
        print(self.read_mode_Raman)
        if self.read_mode_Raman == 'Image':
            self.Andor.set_andor_parameter('Image', 1, 1, 1, self.pixel_number, central_row - rows, central_row + rows)
        elif self.read_mode_Raman == 'Single track':
            self.Andor.set_andor_parameter('SingleTrack', central_row, rows*2)
    
    def Save_Raman(self):
        # for saving single measurement
        
#        if self.Filename_text.toPlainText() != 'File name ...':
#            filename = self.Filename_text.toPlainText()
#        else:
#            filename = 'Andor_data'
        if self.Filename_text_Raman.text().strip() != 'File name ...':
            filename = self.Filename_text_Raman.text().strip()
        else:
            filename = 'Andor_data'

        activeRamandata = self.Datagroup_Raman.create_dataset(filename, data = spectrum_Raman)

        activeRamandata.attrs.create('Description', self.Description_text_Raman.toPlainText())        
        activeRamandata.attrs.create('Exposure_time', round(self.Andor.Exposure, 2))
        activeRamandata.attrs.create('wavelengths', wavelength_Raman)
        activeRamandata.attrs.create('background', 0)
        
        """ save data into text file"""
        now = datetime.now()
        current_date = now.strftime("%y%m%d")
        
        self.folder2save = os.path.join('C:\\Users\\Lab Di Martino\\Documents\\data\\' + self.Andor.CRSID)
        
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
        saveData = np.vstack(np.transpose((\
                                           np.append([0], wavelength_Raman), \
                                           np.append([round(self.Andor.Exposure, 2)], spectrum_Raman))))
        np.savetxt(afile, saveData)
        afile.close()


        """ END: save data into text file"""
       
    def update_figure_lim_Raman(self):
        
        self.Andor.laser_wl = float(self.LaserWavelength_value_Raman.value())
        self.Spectra_display_Raman.fontsize = 16
        def WL2WN(wavelength):
            return 1/(self.Andor.laser_wl*1e-7) - 1/(wavelength*1e-7)
        
        def WN2WL(k):
            return 1e7/(1/(self.Andor.laser_wl*1e-7) - k)
        
        self.Spectra_display_Raman.canvas.ax.set_ylabel('Intensity (counts)', fontsize = self.Spectra_display_Raman.fontsize, color = 'white')
        self.Spectra_display_Raman.canvas.ax.set_xlabel('Wavelength (nm)', fontsize = self.Spectra_display_Raman.fontsize, color = 'white')
        self.Spectra_display_Raman.canvas.ax.tick_params(axis = "y", direction = "in")
        self.Spectra_display_Raman.canvas.ax.tick_params(axis = "x", direction = "in", top = False)
        self.Spectra_display_Raman.canvas.ax.set_xlim([min(wavelength_Raman), max(wavelength_Raman)])
        
        self.Spectra_display_Raman.minY = min(float(self.YLim1_box_Raman.value()), float(self.YLim2_box_Raman.value()))
        self.Spectra_display_Raman.maxY = max(float(self.YLim1_box_Raman.value()), float(self.YLim2_box_Raman.value()))
        if self.Spectra_display_Raman.minY < self.Spectra_display_Raman.maxY:
            self.Spectra_display_Raman.canvas.ax.set_ylim([self.Spectra_display_Raman.minY, self.Spectra_display_Raman.maxY])
        
        """ Display WN axis"""
        
        self.Spectra_display_Raman.minX = min(float(self.XLim1_box_Raman.value()), float(self.XLim2_box_Raman.value()))
        self.Spectra_display_Raman.maxX = max(float(self.XLim1_box_Raman.value()), float(self.XLim2_box_Raman.value()))
        if self.Spectra_display_Raman.minX == self.Spectra_display_Raman.maxX:
            self.Spectra_display_Raman.minX = min(wavelength_Raman)
            self.Spectra_display_Raman.maxX = max(wavelength_Raman)
            
        self.Spectra_display_Raman.canvas.ax.set_xlim([self.Spectra_display_Raman.minX, self.Spectra_display_Raman.maxX])
        
        self.Spectra_display_Raman.canvas.ax2.set_xlabel('Wavenumber (cm$^{-1}$)', 
                                                         fontsize = self.Spectra_display_Raman.fontsize, color = 'white')
#        self.Spectra_display.canvas.ax2.set_xlim([self.Spectra_display.minX, self.Spectra_display.maxX])
        
        kmin = WL2WN(self.Spectra_display_Raman.minX)
        kmax = WL2WN(self.Spectra_display_Raman.maxX)
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
        self.Spectra_display_Raman.canvas.ax2.set_xticks(kticksInNm)
        self.Spectra_display_Raman.canvas.ax2.set_xticklabels(kticksstr)
        self.Spectra_display_Raman.canvas.ax2.set_xlim([self.Spectra_display_Raman.minX, self.Spectra_display_Raman.maxX])
        """ Display WN axis"""
        
        for label in (self.Spectra_display_Raman.canvas.ax.get_xticklabels() +  \
                      self.Spectra_display_Raman.canvas.ax.get_yticklabels() + \
                      self.Spectra_display_Raman.canvas.ax2.get_xticklabels()):
            label.set_fontsize(self.Spectra_display_Raman.fontsize)
        self.Spectra_display_Raman.canvas.fig.tight_layout()
        self.Spectra_display_Raman.canvas.draw()
        
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
        self.OOspectrometer = OOspectrometer
    
    def run(self):
        #TEMPORARY PATCH 26TH APRIL - READ SPECTRUM TWICE, SO THAT ANY CHANGES JUST BEFORE 'TAKE SPECTRUM' WAS CLICKED DONT EFFECT SPECTRUM COLLECTED
        self.current_spec = self.OOspectrometer.read_spectrum()
        self.current_spec = self.OOspectrometer.read_spectrum()
        
        self.OOspectrometer.OO_occupied = 0
        print('Finished taking spec, using int t = '+str(self.OOspectrometer.integration_time))
    
"""start: helper classes for SMC100 stage"""
class sub_SMC100_move(QRunnable):
    def __init__(self, SMC100, axis, value):
        super().__init__()
        self.SMC100 = SMC100
        self.axis = axis
        self.value = value        
        
    def run(self):
        print(' --> Moving SMC100')
        #execute movement
        initial_position = self.SMC100.get_position(self.axis)
        target_position= float(initial_position[0]) + self.value
        self.SMC100.move(target_position, self.axis, relative = False)
        final_position = self.SMC100.get_position(self.axis)
        #mark stage as no longer in use
        global stage_in_use
        stage_in_use = False 
        #print updated info
        print('Position axis ' + str(self.axis) + ': ' + str(final_position) + ' mm')
        
        
class sub_SMC100_move_mid(QRunnable):
    def __init__(self, SMC100):
        super().__init__()
        self.SMC100 = SMC100

    def run(self):
        print('Moving to midpoint of stage range of motion.')
        self.SMC100.move(6, '1', relative=False)
        self.SMC100.move(6, '2', relative=False)
        self.SMC100.move(6, '3', relative=False)
        global stage_in_use
        stage_in_use = False 
        print('All axes of stage moved to position = 6 mm.')
        print(' ')
        
        
class sub_initialise_stage(QRunnable):
    def __init__(self, SMC100):
        super().__init__()
        self.SMC100 = SMC100

    def run(self):
        print(' ')
        print('Wait for stage initialisation...')
#        self.SMC100.reset_and_configure()
        self.SMC100.home()
        global stage_in_use
        stage_in_use = False
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
    def __init__(self, OOspectrometer, Smaract_stage, position, channel):
        super().__init__()
        # this is local access to Smaract_stage. "self" here is the "sub_stage_on_hold"
        self.Smaract_stage = Smaract_stage
        self.channel = channel
        self.position = position
        self.OOspectrometer = OOspectrometer
    
    def run(self):
        
        while self.OOspectrometer.OO_occupied > 0:
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
        self.MiTC_handle.setSampleTemp('Sample', self.setSampleTemp)
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
         if str(self.MiPS_handle.getMiPSHeater) == 'OFF':
             self.MiPS_handle.setMiPSHeater('Magnet', 'ON')
             # free the communication
             MiPS_communication = 0
             self.MiPSHeater_button.setText('Heater ... on')
             time.sleep(30)
             self.MiPSHeater_button.setText('Heater ON')
         elif str(self.MiPS_handle.getMiPSHeater) == 'ON':
             self.MiPS_handle.setMiPSHeater('Magnet', 'OFF')
             # free the communication
             MiPS_communication = 0
             self.MiPSHeater_button.setText('Heater ON')
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

""" Add parts for controlling the Raman spectrometer Trung 12.04.24 ----------------------"""
""" Name of all boxes, plot from here are local parameters, so don't need to have
the suffix _Raman. But they are the items with suffix _Raman in the UI """ 
""" Only the spectrum_Raman and wavelength_Raman needs to be correct since they
are global parameters """
class sub_get_temperature_Raman(QRunnable):
    def __init__(self, Andor_temp_button, Andor):
        super().__init__()
        self.Andor_temp_button = Andor_temp_button
        self.Andor = Andor
        
    def run(self):
        current_Andor_temp = float(self.Andor.CurrentTemperature)
        self.Andor_temp_button.setText('Andor temp. ' + str(current_Andor_temp) + ' C')
        
        self.Andor.Andor_occupied = 0
        
class sub_capture_Raman(QRunnable):
    def __init__(self, Andor, Spectra_display, Camerastate_output, YLim1_box, YLim2_box):
        super().__init__()
        self.Andor = Andor
        self.Spectra_display = Spectra_display
        self.Camerastate_output = Camerastate_output
        self.laser_wl = self.Andor.laser_wl
        self.YLim1_box = YLim1_box
        self.YLim2_box = YLim2_box
        
    def run(self):
        global spectrum_Raman
        global wavelength_Raman
        
        self.Camerastate_output.setStyleSheet('background-color: rgb(255, 140, 140);')
        
        imageArray, num_of_images, image_shape = self.Andor.capture()
        imageArray.reverse()
        
        spectrum_Raman = np.array(imageArray)
        
        self.Spectra_display.canvas.ax.cla()

        def WL2WN(wavelength):
            return 1/(self.laser_wl*1e-7) - 1/(wavelength*1e-7)
        
        def WN2WL(k):
            return 1e7/(1/(self.laser_wl*1e-7) - k)

        p = self.Spectra_display.canvas.ax.plot(wavelength_Raman, spectrum_Raman)
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

        self.Andor.Andor_occupied = 0
        
        
class sub_capture_and_save_Raman(QRunnable):
    def __init__(self, Datagroup, Andor, Spectra_display, Camerastate_output, \
                 Repeat_box_value, Filename_text, Description_text, YLim1_box, YLim2_box, path2save):
        
        super().__init__()
        self.Datagroup = Datagroup
        self.Andor = Andor
        self.Spectra_display = Spectra_display
        self.Camerastate_output = Camerastate_output
        self.laser_wl = self.Andor.laser_wl
        self.Repeat_box_value = Repeat_box_value
        self.Filename_text = Filename_text
        self.Description_text = Description_text
        self.YLim1_box = YLim1_box
        self.YLim2_box = YLim2_box
        self.path2save = path2save
        
    def run(self):
        global spectrum_Raman
        global wavelength_Raman
        
        i = 0
        while i <= self.Repeat_box_value:
            print(i)
            self.Camerastate_output.setText('Capturing ... ' + str(i))
            self.Camerastate_output.setStyleSheet('background-color: rgb(255, 140, 140);')
            imageArray, num_of_images, image_shape = self.Andor.capture()
            imageArray.reverse()
            
            spectrum_Raman = np.array(imageArray)
           
            self.Spectra_display.canvas.ax.cla()

            
            wavelength2 = np.array(wavelength_Raman)
            def WL2WN(wavelength):
                return 1/(self.laser_wl*1e-7) - 1/(wavelength*1e-7)
            
            def WN2WL(k):
                return 1e7/(1/self.laser_wl*1e-7 - k)
            
            def WN2WL(k):
                return 1e7/(1/(self.laser_wl*1e-7) - k)

            k = WL2WN(wavelength2)
    
                
            """ save file """
#            if self.Filename_text.toPlainText() != 'File name ...':
#                filename = self.Filename_text.toPlainText()
#            else:
#                filename = 'Andor_data'
            if self.Filename_text.text().strip() != 'File name ...':
                filename = self.Filename_text.text().strip()
            else:
                filename = 'Andor_data'
            activeRamandata = self.Datagroup.create_dataset(filename, data = spectrum_Raman)
    
            activeRamandata.attrs.create("Description", self.Description_text.toPlainText())        
            activeRamandata.attrs.create("Exposure_time", round(self.Andor.Exposure, 2))
            activeRamandata.attrs.create("wavelengths", wavelength_Raman)
            activeRamandata.attrs.create("background", 0)
            i += 1
            
            # Save as increment file name
            increment_ = 0
            
            fullfile = self.path2save + '\\Raman_' + filename
            if os.path.exists(fullfile + '.txt'):
                while os.path.exists(fullfile + '_' + str(increment_) + '.txt'):
                    increment_ = increment_ + 1
                print('Save file nr. ' + str(increment_))
                outputfile = fullfile + '_' + str(increment_) + '.txt'
            else:
                outputfile = fullfile + '.txt'
                
            afile = open(outputfile, 'w')
            saveData = np.vstack(np.transpose((\
                                           np.append([0], wavelength_Raman), \
                                           np.append([round(self.Andor.Exposure, 2)], spectrum_Raman))))
            np.savetxt(afile, saveData)
            afile.close()
            
            
            ''' Plot '''
            p = self.Spectra_display.canvas.ax.plot(wavelength_Raman, spectrum_Raman)
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
            
            
        self.Camerastate_output.setText('Done!')
        self.Camerastate_output.setStyleSheet('background-color: rgb(180, 227, 255);')
        
        self.Andor.Andor_occupied = 0




class sub_Live_mode_Raman(QRunnable):
    def __init__(self, Andor, Spectra_display, Camerastate_output, \
                 YLim1_box, YLim2_box):
        
        super().__init__()

        self.Andor = Andor
        self.Spectra_display = Spectra_display
        self.Camerastate_output = Camerastate_output
        self.laser_wl = self.Andor.laser_wl
        self.YLim1_box = YLim1_box
        self.YLim2_box = YLim2_box
        
    def run(self):
        global spectrum_Raman
        global wavelength_Raman
        global Abort_live_mode
        

        while Abort_live_mode < 1:

            self.Camerastate_output.setText('Capturing ... ')
            self.Camerastate_output.setStyleSheet('background-color: rgb(255, 140, 140);')
            imageArray, num_of_images, image_shape = self.Andor.capture()
            imageArray.reverse()
            
            spectrum_Raman = np.array(imageArray)
           
            self.Spectra_display.canvas.ax.cla()
            
            wavelength2 = np.array(wavelength_Raman)
            def WL2WN(wavelength):
                return 1/(self.laser_wl*1e-7) - 1/(wavelength*1e-7)
            
            def WN2WL(k):
                return 1e7/(1/(self.laser_wl*1e-7) - k)

            k = WL2WN(wavelength2)
    
            p = self.Spectra_display.canvas.ax.plot(wavelength_Raman, spectrum_Raman)
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
        
        self.Andor.Andor_occupied = 0

class sub_captureImage_Raman(QRunnable):
    def __init__(self, Andor, Spectra_display, Camerastate_output, YLim1_box, YLim2_box):
        super().__init__()
        self.Andor = Andor
        self.Spectra_display = Spectra_display
        self.Camerastate_output = Camerastate_output
        self.laser_wl = self.Andor.laser_wl
        self.YLim1_box = YLim1_box
        self.YLim2_box = YLim2_box
        
    def run(self):
        global spectrum_Raman
        global wavelength_Raman
        
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
        real_x = wavelength_Raman
        
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

        
        self.Andor.Andor_occupied = 0
        



































