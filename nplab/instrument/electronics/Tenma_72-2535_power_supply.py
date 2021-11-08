# -*- coding: utf-8 -*-
"""
Created on Mon Nov  8 15:00:00 2021

@author: dk515
"""

from builtins import str
from nplab.instrument.serial_instrument import SerialInstrument
import serial
import time

class Tenma_72_2535(SerialInstrument):
    """
    Class to control a Tenma 72-2535 programmable power supply. Setting voltage,
    current and output works fine, but other commands sometimes time out or crash,
    so need further testing.
    """
    
    port_settings = dict(baudrate=9600,
                         bytesize=serial.EIGHTBITS,
                         parity=serial.PARITY_NONE,
                         stopbits=serial.STOPBITS_ONE,
                         #xonxoff=True,
                         timeout= 1.0)
    termination_character = '\n'
    
    def __init__(self, port = None):
        SerialInstrument.__init__(self, port=port)
        
    def IDN(self):
        self.write('*IDN?')
        time.sleep(0.2)
        return self.readline()
    
    def set_voltage(self, voltage):
        """ Set voltage in V """
        self.write('VSET1:' + str(voltage))
        
    def set_current(self, current):
        """ Set current in A """
        self.write('ISET1:' + str(current))
        
    def read_current(self):
        """ Returns actual output current """
        self.write('IOUT1?')
        time.sleep(0.5)
        return self.readline()
    
    def read_voltage(self):
        """ Returns actual output voltage """
        self.write('VOUT1?')
        time.sleep(0.5)
        return self.readline()
    
    def output(self, output):
        """ Turn output on/off, 0: OFF, 1: ON """
        self.write('OUT' + str(output))
      
if __name__ == "__main__":
    my_Tenma = Tenma_72_2535(port = 'COM4')