# -*- coding: utf-8 -*-
"""
Created on Sat Jul 08 19:47:22 2017

@author: Hera
"""
from nplab.instrument.camera.Andor import Andor, AndorUI
from nplab.instrument.spectrometer.Kymera import Kymera
from nplab.utils.notified_property import NotifiedProperty
#import numpy as np
class Kandor(Andor):
    ''' Wrapper class for the kymera and the andor
    '''
    
    def __init__(self, pixel_number=1024,
                 pixel_width=26,
                 use_shifts=False, 
                 laser_wl=632.8,
                 
                 white_shutter=None):
        
        super().__init__()
        self.kymera = Kymera()
        self.kymera.pixel_number = pixel_number
        self.kymera.pixel_width = pixel_width
        self.use_shifts = use_shifts
        self.laser_wl = laser_wl
        self.white_shutter = white_shutter
        self.kymera.center_wavelength = 633.0
        self.metadata_property_names += ('slit_width', 'wavelengths')
        self.ImageFlip = 0
        

    def get_x_axis(self):
        return self.kymera.GetCalibration()[::-1]
    x_axis = property(get_x_axis) #This is grabbed by the Andor code 

    @property
    def slit_width(self):
        return self.kymera.slit_width
    
    @property 
    def wavelengths(self):
        return self.get_x_axis()
    
    def Capture(_AndorUI):
        _AndorUI.Andor.raw_image(update_latest_frame = True)
    setattr(AndorUI, 'Capture', Capture)

if __name__ == '__main__':
    k = Kandor()
    k.show_gui(block = False)
    ky = k.kymera
    ky.show_gui(block = False)
#    k.MultiTrack = (2, 3, 50)