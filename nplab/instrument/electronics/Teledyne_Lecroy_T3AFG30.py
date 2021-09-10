# -*- coding: utf-8 -*-
"""
Created on Mon Sep  6 13:59:17 2021

@author: dk515, aj619, trr30
"""

"""
Class for Teledyne Lecroy T3AFG30 arbitrary function generator.
"""

from nplab.instrument.visa_instrument import VisaInstrument

class Teledyne_Lecroy_T3AFG30(VisaInstrument):
    """Interface to the T3AFG30 arbitrary function generator."""
    def __init__(self, address='USB0::0xF4EC::0xEE38::T0102C21150130::INSTR'):
        super(Teledyne_Lecroy_T3AFG30, self).__init__(address)
        self.instr.read_termination = '\n'
        self.instr.write_termination = '\n'
        self.reset()
        self.outputs_off()

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
        
    def set_voltage(self, channel, voltage):
        """
        Set amplitude of specific channel in Vpp.
        <channel> 'C1' or 'C2'
        <voltage> '3' for 3 Vpp, e.g.
        """
        self.write(channel + ':BSWV AMP,' + voltage)
        
    def set_frequency(self, channel, freq):
        """
        Set frequency of specific channel in Hz.
        <channel> 'C1' or 'C2'
        <freq> '1000' for 1 kHz, e.g.
        """
        self.write(channel + ':BSWV FRQ,' + freq)
        
    def set_arb_func(self, channel, func_name):
        """
        Assign stored waveform to channel.
        <channel> 'C1' or 'C2'
        <func_name> 'PUND', e.g.
        """
        self.write(channel + ':ARWV NAME,' + func_name)
        'must match exactly, does not throw error if no stored waveform matches'
        
    def set_sync_signal(self, channel, state, src):
        """
        Turns on external sync signal based on channel or its modulation function
        <channel> 'C1' or 'C2'
        <state> 'ON' or 'OFF'
        <src> 'CH1' 'CH2' 'MOD_CH1' or 'MOD_CH2'
        """
        self.write(channel + ':SYNC ' + state + ',TYPE,' + src)
        
    def outputs_off(self):
        """
        Turns of CH1 and CH2 output
        """
        self.set_output('C1','OFF')
        self.set_output('C2','OFF')
     
if __name__ == '__main__':
    myT3AFG30 = Teledyne_Lecroy_T3AFG30()