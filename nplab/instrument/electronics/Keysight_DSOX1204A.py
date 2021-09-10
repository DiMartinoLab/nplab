# -*- coding: utf-8 -*-
"""
Created on Mon Sep  6 15:13:41 2021

@author: dk515, aj619, trr30
"""

"""
Class for Keysight DSOX1204A digital oscilloscope.
"""

from nplab.instrument.visa_instrument import VisaInstrument
import sys
from comtypes.client import GetModule
from comtypes.client import CreateObject
#if not hasattr(sys, "frozen"):
#    GetModule("C:\Program Files (x86)\IVI Foundation\VISA\VisaCom\GlobMgr.dll")
import comtypes.gen.VisaComLib as VisaComLib
#import numpy as np


class Keysight_DSOX1204A(VisaInstrument):
    """Interface to the DSOX1204A oscilloscope."""
    def __init__(self, address='USB0::0x2A8D::0x0386::CN60476268::0::INSTR'):
        super(Keysight_DSOX1204A, self).__init__(address)
        self.instr.read_termination = '\n'
        self.instr.write_termination = '\n'
        self.clear()
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
            
            
    def trigger(self,trigger_pattern):
        self.write(':TRIGger:PATTern {}'.format(trigger_pattern)) # command is :TRIGger:PATTern <pattern>
            
    def autoscale(self):
        self.write(':AUToscale') #The :AUToscale command evaluates all input signals and sets the correct
                                #conditions to display the signals
            
    def acquire_type(self,control_value):
        self.write(':ACQuire:TYPE {}'.format(control_value)) #select normal, avergae or high resolution data
        
    def acquire_mode(self,set_mode):
        self.write(':ACQuire:MODE {}'.format(set_mode)) # other options are segmented. Option selected is RealTime
        
        
#The :DIGitize command is a specialized RUN command. It causes the instrument
#to acquire waveforms according to the settings of the :ACQuire commands
#subsystem. When the acquisition is complete, the instrument is stopped.        
    def digitize(self,channel_numbers_on):
        for items in channel_numbers_on:
            self.write(':DIGitize CHANnel{}'.format(items))
            
#            
#    def single_run(self):
#        self.write(':SINGle')  #The :SINGle command causes the instrument to acquire a single trigger of data.
        
    def run(self):
        self.write(':RUN') #The :RUN command starts repetitive acquisitions
        
#    def stop(self):
#        self.write(':STOP') #The :STOP command stops the acquisition.
        
    def do_query_string(query):
        myScope.WriteString("%s" % query, True)
        result = myScope.ReadString()
       # check_instrument_errors(query)
        return result
        
        
    def waveform_source(self):
        self.write(':WAVeform:SOURce CHANnel1')  
        


    def download_waveform(self):
        self.write(':WAVeform:POINts:MODE RAW')

          
        
if __name__ == '__main__':
    
    rm = CreateObject("VISA.GlobalRM", \
                      interface=VisaComLib.IResourceManager)
    myScope = CreateObject("VISA.BasicFormattedIO", \
                           interface=VisaComLib.IFormattedIO488)
    myScope.IO = \
    rm.Open("TCPIP0::141.121.237.208::hislip0::INSTR")
    
    
    myDSOX1204A = Keysight_DSOX1204A()
    myDSOX1204A.channel_display_on(['1','2'])
    myDSOX1204A.channel_display_off(['3','4'])
    myDSOX1204A.acquire_type('NORMal') # other options are 'AVERage','HRESolution','PEAK'
    myDSOX1204A.acquire_mode('SEGMented') # other options are segmented. Option selected is RealTime
    myDSOX1204A.digitize(['1'])
    myDSOX1204A.save_waveform()
   # myDSOX1204A.trigger('\"1\"[,CHANnel1,POSitive]')
    myDSOX1204A.autoscale()
    myDSOX1204A.download_waveform('50')
    myDSOX1204A.download_waveform()
    myDSOX1204A.do_query_string()
    
    qresult = myDSOX1204A.do_query_string(":WAVeform:POINts:MODE?")
    print("Waveform points mode: %s" % qresult)
    
    