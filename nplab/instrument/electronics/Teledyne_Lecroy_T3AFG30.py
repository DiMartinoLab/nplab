# -*- coding: utf-8 -*-
"""
Created on Mon Sep  6 13:59:17 2021

@author: dk515, aj619, trr30
"""

"""
Class for Teledyne Lecroy T3AFG30 arbitrary function generator.
"""

from nplab.instrument.visa_instrument import VisaInstrument
import binascii
import time as time
import numpy as np


class Teledyne_Lecroy_T3AFG30(VisaInstrument):
    """Interface to the T3AFG30 arbitrary function generator."""
    def __init__(self, address='USB0::0xF4EC::0xEE38::T0102C21150130::INSTR'):
        super(Teledyne_Lecroy_T3AFG30, self).__init__(address)
        self.instr.read_termination = '\n'
        self.instr.write_termination = '\n'
        self.reset()

    def reset(self):
        """Reset the instrument to its default state."""
        self.write('*RST')
        
    def set_output(self, channel, state):
        """
        Set ouput state of channels.
        :channel: 'C1' or 'C2'
        :state: 'ON' or 'OFF'
        """
        self.write(channel + ':OUTP ' + state)
        
    def create_wave_file(self, wave_points):
        """create a file"""
        f = open("wave1.bin", "wb")
        for a in wave_points:
            b = hex(a)
            b = b[2:]
            len_b = len(b)
            if (0 == len_b):
                b = '0000'
            elif (1 == len_b):
                b = '000' + b
            elif (2 == len_b):
                b = '00' + b
            elif (3 == len_b):
                b = '0' + b
            print(b)
            #b = b[2:4] + b[:2]
            #print(b)
            c = binascii.a2b_hex(b) #Hexadecimal integer to ASCii encoded string
            print(c)
            f.write(c)
        f.close()

    def send_wave_data(self):
        """send wave1.bin to the device"""
        f = open("wave1.bin", "rb") #wave1.bin is the waveform to be sent
        data = f.read()
        print(data)
        print('write bytes:',len(data))
        self.write("C1:WVDT WVNM,wave1,FREQ,20.0,AMPL,4.0,OFST,0.0,PHASE,0.0,WAVEDATA,%s" %(data))
    #T3AFG5, 10 and T3AFG40, 80, 120
        self.write("C1:ARWV NAME,wave1")
        f.close()
        
    def get_wave_data(self):
        """get wave from the device"""
        f = open("wave2.bin", "wb") #save the waveform as wave2.bin
        self.write("WVDT? user,test") #T3AFG5, 10 and T3AFG40, 80, 120
        time.sleep(1)
        data = self.read()
        print(data)
        data_pos = data.find("WAVEDATA,") + len("WAVEDATA,")
        print(data_pos)
        print(data[0:data_pos])
        wave_data = data[data_pos:]
        print('read bytes:',len(wave_data))
        f.write(wave_data)
        f.close()
    
    
            
if __name__ == '__main__':
    myT3AFG30 = Teledyne_Lecroy_T3AFG30()
    print('test2')
    
    wave_points = [0x0010, 0x0020, 0x0030, 0x0040, 0x0050, 0x0060, 0x0070, 0xff7f]
    #myT3AFG30.create_wave_file(wave_points)
    d = ''
    e = ''
    
    for a in wave_points:
            b = hex(a)
            b = b[2:]
            len_b = len(b)
            if (0 == len_b):
                b = '0000'
            elif (1 == len_b):
                b = '000' + b
            elif (2 == len_b):
                b = '00' + b
            elif (3 == len_b):
                b = '0' + b
            d = np.append(d, b)
            print(d)
            e = np.append(e, binascii.a2b_hex(b))
            print(e)