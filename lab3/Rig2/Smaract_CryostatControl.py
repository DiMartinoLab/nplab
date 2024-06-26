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

class Smaract_CryostatControl(Experiment, QtWidgets.QWidget, UiTools):
    #inherit QWidget methods
    def __init__(self, Stage_instance, \
                 threadCount, pool, \
                 ui_file = os.path.join(os.path.dirname(__file__),'Smaract_CryostatControl.ui'),   parent=None):
        
        #check wch stage is connected

        self.Smaract_stage = Stage_instance
        
        # make this window available for the main program
        super(Smaract_CryostatControl, self).__init__() 
        uic.loadUi(ui_file, self)
        


        """setup stage controls and display - using olympus OR smaract"""
        self.disconnect_stage_button.clicked.connect(self.disconnect_stage)
        
        # dmk50 April 2023
        
        #buttons and display for z setpoint
        self.Set_New_Coord_origin_button.clicked.connect(self.Set_New_Origin)
        self.Go_To_Stage_Origin_button.clicked.connect(self.Go_To_Stage_Origin)
        self.origin_stored_checkbox.setChecked(False)
        
            
            
        
        # Trung Nov. 2022
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
    

    
        """ parallel processing"""
#        self.threadCount = QThreadPool.globalInstance().maxThreadCount()
#        self.pool = QThreadPool.globalInstance()
        self.threadCount = threadCount
        self.pool = pool

        
    # Trung Nov. 2022
    # added parts for the SmarAct stage control
    # Further infos and code: C:\SmarAct\MCS2\SDK\Python\examples
    # self.Smaract_stage = d_handle
  
    
    """START: methods for control of stages"""
    def disconnect_stage(self):
        smaract_package.Close(self.Smaract_stage)
        print('Smaract stage disconnected.')
        
            
    def Set_New_Origin(self):
        #this will take current coordinates of the three smaracts and save them as 'origin' - note that correct frame of reference must be used

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
        


        
    def Go_To_Stage_Origin(self):
        
            
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

    def print_stage_positions(self):

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
        
   

    def make_window(self):
        app = get_qt_app()
        self.show()
        app.exec_()
        return self


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
        self.OOspectrometer = OOspectrometer
    
    def run(self):
        
        while self.OOspectrometer.OO_occupied > 0:
            cposition = smaract_package.GetProperty_i64(self.Smaract_stage, self.channel, smaract_package.Property.POSITION)
            if np.abs(cposition - self.position) > 500000: # = 500nm
                smaract_package.Move(self.Smaract_stage, self.channel, self.position, 0)
                print('Stage position corrected')
"""end: helper class for smaract stage"""  


