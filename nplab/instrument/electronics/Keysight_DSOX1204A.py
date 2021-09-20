# -*- coding: utf-8 -*-
"""
Created on Mon Sep  6 15:13:41 2021

@author: dk515, aj619, trr30
"""

"""
Class for Keysight DSOX1204A digital oscilloscope.
"""

"""import modules"""


from nplab.instrument.visa_instrument import VisaInstrument
import pyvisa
import struct
#import matplotlib.pyplot as plt
import numpy as np
import time

class Keysight_DSOX1204A(VisaInstrument):
    """Interface to the DSOX1204A oscilloscope."""
    def __init__(self, address='USB0::0x2A8D::0x0386::CN60476268::0::INSTR'):
        super(Keysight_DSOX1204A, self).__init__(address)
        self.instr.read_termination = '\n'
        self.instr.write_termination = '\n'
       # self.clear()
       # self.reset()
       # self.data_format('ASCII')
       # self.set_waveform_points('2000')
       # self.channel_attenuation_1x('1')
       # self.channel_attenuation_1x('2')
       # self.channel_attenuation_1x('3')
       # self.channel_attenuation_1x('4')
        
    def reset(self):
        """Reset the instrument to its default state."""
        self.write('*RST')
        
    def clear(self):
        """The *CLS common command clears the status data structures"""
        self.write('*CLS')
        
    def channel_display_on(self,channels):
        """Turns on channels for display/capturing.
        <channels> ['1'] or ['1','2'], i.e.
        """
        for item in channels:
            self.write(':CHANnel{}:DISPlay 1'.format(item))
            
    def channel_display_off(self,channels):
        """Turns off channels for display/capturing.
        <channels> ['1'] or ['1','2'], i.e.
        """
        for item in channels:
            self.write(':CHANnel{}:DISPlay 0'.format(item))
    
    def channel_range(self, channel, vertical_range):
        self.write('CHANnel' + channel + ':RANGe ' + vertical_range)   
        
    def channel_offset(self, channel, offset):
        self.write('CHANnel' + channel + ':OFFSet ' + offset)
        
    def channel_attenuation_1x(self, channel):
        self.write('CHANnel' + channel + ':PROBe 1')
        
    def channel_invert_toggle(self, channel):
        #Note: this is for visual purposes only! Does not apply to raw values.
        state = int(self.query(':CHANnel' + channel + ':INVert?'))
        self.write(':CHANnel' + channel + ':INVert ' + str((state + 1) % 2))
        
    def channel_invert(self, channel, state):
        #Note: this is for visual purposes only! Does not apply to raw values.
        """
        <state> 'ON' or 'OFF; '1' or '0'
        """
        self.write(':CHANnel' + channel + ':INVert ' + state)
            
    def trigger_edge_source(self, source):
        """Sets source of trigger signal.
        <source> 'CHAN<n>', 'EXT', 'LINE', 'WGEN', 'NONE'
        <n> for channels 1-4
        """
        self.write(':TRIGger:EDGE:SOURce ' + source) # command is :TRIGger:PATTern <pattern>    
        
    def trigger_edge(self):
        self.write(':TRIGger:MODE EDGE')
        
    def trigger_edge_polarity(self):    
        self.write(":TRIGger:EDGE:SLOPe POSitive")
    
    def trigger_sweep(self, sweep):
        """
        Sets trigger sweep--AUTO refreshes without trigger, NORMal waits for new trigger.
        <sweep> 'AUTO' or 'NORMal'
        """
        self.write(':TRIGger:SWEep ' + sweep)
        
    def trigger_level_external(self, level):
        self.write(':EXTernal:LEVel ' + level)
        
    def timebase_reference(self, reference):
        """
        Aligns time reference to 1 div left/right or in center of screen.
        <reference> 'LEFT', 'CENTer', 'RIGHt'
        """
        self.write(':TIMebase:REFerence ' + reference)
        
    def timebase_range(self, rang):
        """
        Sets the full-scale horizontal time in seconds for the main window.
        The range is 10 times the current time-per-division setting.
        <rang> '<time for 10 div in sec>', i.e., '1e-3'
        """
        self.write(':TIMebase:RANGe ' + rang)
        
    def timebase_position(self, position):
        """
        Sets the time interval between the trigger event and the display reference point on the screen.
        <position> '<NR3 format number>', i.e., '100E-06'
        """
        self.write(':TIMebase:POSition ' + position)
        
        #TODO
        #set external trigger threshold voltage (2V)
        
    def autoscale(self):
        self.write(':AUToscale') #The :AUToscale command evaluates all input signals and sets the correct
                                #conditions to display the signals
                                
    def autoscale_custom(self, channel, frequency, amplitude):
        self.timebase_reference('LEFT')
        self.timebase_range(str(1/frequency))
        self.timebase_position(str(0.1/frequency))
        self.channel_range(str(channel),str(amplitude))
        
            
    def acquire_type(self,control_value):
        self.write(':ACQuire:TYPE {}'.format(control_value)) #select normal, avergae or high resolution data
        
    def acquire_mode(self,set_mode):
        self.write(':ACQuire:MODE {}'.format(set_mode)) # other options are segmented. Option selected is RealTime
        
                
    def vertical_scale(self,vertical_scale):
        self.write(':CHANnel1:SCALe {}'.format(vertical_scale))
        
    def vertical_range(self,vertical_range):
        self.write(':CHANnel1:RANGe ' + vertical_range)        
    
    def digitize(self,channel_numbers_on):  #specialized RUN command
        for items in channel_numbers_on:
            self.write(':DIGitize CHANnel{}'.format(items))
            
            
    def single_run(self):
        self.write(':SINGle')  #The :SINGle command causes the instrument to acquire a single trigger of data.
        
    
    def run(self):
        self.write(':RUN') #The :RUN command starts repetitive acquisitions
        
    def stop(self):
        self.write(':STOP') #The :STOP command stops the acquisition.
        
        
        
    def waveform_mode(self, mode): # sets waveform mode
        """
        <mode> 'NORMal', 'MAXimum', 'RAW'
        """
        self.write(':WAVeform:POINts:MODE ' + mode)
    
    def set_waveform_points(self, num_points): # sets number of points to acquire
        self.write(':WAVeform:POINts ' + num_points)
        
    def waveform_source(self,waveform_channel_source):
        self.write(':WAVeform:SOURce CHANnel{}'.format(waveform_channel_source))
        
        
    def waveform_data_format(self, data_format):
        self.write(':WAVEFORM:FORMAT ' + data_format)

    
if __name__ == '__main__':
    
    rm = pyvisa.ResourceManager() # call resource manager from module pyvisa
    rm.list_resources() #lists all resources available in case one forgets the resoucrce address
    myScope = rm.open_resource('USB0::0x2A8D::0x0386::CN60476268::0::INSTR') #opens oscilloscope
    #myScope.timeout = 0 #to set a timeout for the device operation
    myDSOX1204A = Keysight_DSOX1204A()
    
    test = False
    
    if test:
        #myDSOX1204A.autoscale_custom('1',freq,amp)
        myT3AFG30.set_output('C1','ON')
        time.sleep(0.150)
        myDSOX1204A.single_run()
        time.sleep(0.1)
        #time.sleep(0.2)
        myT3AFG30.set_output('C1','OFF')
        #time.sleep(1)
        myDSOX1204A.waveform_source('1')
        c1_data = myDSOX1204A.query(':WAV:DATA?')
        myDSOX1204A.waveform_source('2')
        c2_data = myDSOX1204A.query(':WAV:DATA?')
    
        x_increment = float(myDSOX1204A.query(":WAVeform:XINCrement?"))
        x_origin = float(myDSOX1204A.query(":WAVeform:XORigin?"))
        y_increment = float(myDSOX1204A.query(":WAVeform:YINCrement?"))
        y_origin = float(myDSOX1204A.query(":WAVeform:YORigin?"))
        y_reference = float(myDSOX1204A.query(":WAVeform:YREFerence?"))
    
        x_split = c1_data.split(',')
        x2_split = c2_data.split(',')
    
        f = open("waveform_data.csv", "w")
        for i in range(0, len(x_split) - 1):
            if i is 0:
                x_sep = x_split[i][11:]
                x2_sep = x2_split[i][11:]
            else:
                x_sep = x_split[i][1:]
                x2_sep = x2_split[i][1:]
            time_val = x_origin + (i * x_increment)
            #voltage = (((float(x_sep) - y_reference) * y_increment) + y_origin)
            voltage1 = float(x_sep)
            voltage2 = float(x2_sep)
            f.write("%E, %f, %f\n" % (time_val, voltage1, voltage2))
            #f.write(x_sep + '\n')
        f.close()


    



    

    

    

    
    