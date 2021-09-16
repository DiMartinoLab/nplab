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
import matplotlib.pyplot as plt
import numpy as np
  



class Keysight_DSOX1204A(VisaInstrument):
    """Interface to the DSOX1204A oscilloscope."""
    def __init__(self, address='USB0::0x2A8D::0x0386::CN60476268::0::INSTR'):
        super(Keysight_DSOX1204A, self).__init__(address)
        self.instr.read_termination = '\n'
        self.instr.write_termination = '\n'
       # self.clear()
        self.reset() 
        
    def reset(self):
        """Reset the instrument to its default state."""
        self.write('*RST')
        
    def clear(self):
        """The *CLS common command clears the status data structures"""
        self.write('*CLS')
        
    def channel_display_on(self,channel_numbers_on):
        for item in channel_numbers_on:
            self.write(':CHANnel{}:DISPlay 1'.format(item))
            
    def channel_display_off(self,channel_numbers_off):
        for item in channel_numbers_off:
            self.write(':CHANnel{}:DISPlay 0'.format(item))
            
            
    def trigger_source(self):
        self.write(':TRIGger:EDGE:SOURce CHANnel1') # command is :TRIGger:PATTern <pattern>
        
        
    def trigger_edge(self):
        self.write(':TRIGger:MODE EDGE')
        
    def trigger_edge_polarity(self):    
        self.write(":TRIGger:EDGE:SLOPe POSitive")
            
    def autoscale(self):
        self.write(':AUToscale') #The :AUToscale command evaluates all input signals and sets the correct
                                #conditions to display the signals
            
    def acquire_type(self,control_value):
        self.write(':ACQuire:TYPE {}'.format(control_value)) #select normal, avergae or high resolution data
        
    def acquire_mode(self,set_mode):
        self.write(':ACQuire:MODE {}'.format(set_mode)) # other options are segmented. Option selected is RealTime
        
                
    def set_vertical_scale(self,vertical_scale_value):
        self.write(':CHANnel1:SCALe {}'.format(vertical_scale_value))
    
    def digitize(self,channel_numbers_on):  #specialized RUN command
        for items in channel_numbers_on:
            self.write(':DIGitize CHANnel{}'.format(items))
            
            
    def single_run(self):
        self.write(':SINGle')  #The :SINGle command causes the instrument to acquire a single trigger of data.
    
    def run(self):
        self.write(':RUN') #The :RUN command starts repetitive acquisitions
        
    def stop(self):
        self.write(':STOP') #The :STOP command stops the acquisition.
        
        
        
    def waveform_mode(self): # sets waveform mode
        self.write(':WAVeform:POINts:MODE RAW')
    
    def set_waveform_points(self): # sets number of points
        self.write(':WAVeform:POINts 1800')
        
    def waveform_source(self,waveform_channel_source):
        self.write(':WAVeform:SOURce CHANnel{}'.format(waveform_channel_source))
        
        
    def waveform_data_format(self):
        self.write(':WAVEFORM:FORMAT BYTE')

    
if __name__ == '__main__':
    
    rm = pyvisa.ResourceManager() # call resource manager from module pyvisa
    rm.list_resources() #lists all resources available in case one forgets the resoucrce address
    myScope = rm.open_resource('USB0::0x2A8D::0x0386::CN60476268::0::INSTR') #opens oscilloscope
    myScope.timeout = 25000 #to set a timeout for the device operation

     
    myDSOX1204A = Keysight_DSOX1204A()
  # myDSOX1204A.autoscale() # optional, creates a time lag for other instructions to execute and throws error
    myDSOX1204A.channel_display_on(['1']) # channel 1 display alone is on
    myDSOX1204A.trigger_source()
    myDSOX1204A.trigger_edge()
    myDSOX1204A.trigger_edge_polarity()
  #  myDSOX1204A.set_vertical_scale('20')
    myDSOX1204A.acquire_type('NORMal') # other options are 'AVERage','HRESolution','PEAK'
    #myDSOX1204A.run()
    #myDSOX1204A.single_run()
    #myDSOX1204A.digitize(['1']) # digitized only channel 1.
    myDSOX1204A.waveform_mode() # set waveform mode
    myDSOX1204A.set_waveform_points()
    myDSOX1204A.waveform_source(['1'])
    myDSOX1204A.waveform_data_format() # can be set to BYTE, WORD, ASCii
    myDSOX1204A.acquire_mode('RTIMe') # other options are segmented. Option selected is RealTime

    num_cycles = 2
    values_all_cycles = np.array([])
    for cycles in range(num_cycles):    
        values_raw_binary = myScope.query_binary_values('WAV:DATA?', datatype='s')
        values_unpacked = struct.unpack("%dB" % len(values_raw_binary[0]), values_raw_binary[0])
        values_all_cycles = np.append(values_all_cycles,values_unpacked)
        
    x_increment = myDSOX1204A.write(":WAVeform:XINCrement?")[0]
    x_origin = myDSOX1204A.write(":WAVeform:XORigin?")[0]
    y_increment = myDSOX1204A.write(":WAVeform:YINCrement?")[0]
    y_origin = myDSOX1204A.write(":WAVeform:YORigin?")[0]
    y_reference = myDSOX1204A.write(":WAVeform:YREFerence?")[0]
   
    f = open("waveform_data.csv", "w")
    for i in range(0, len(values_all_cycles) - 1):
        time_val = x_origin + (i * x_increment)
        voltage = (((int(values_all_cycles[i]) - y_reference) * y_increment) + y_origin)
        f.write("%E, %f\n" % (time_val, voltage))
    f.close()
    
    plt.plot(values_all_cycles)



    



    

    

    

    
    