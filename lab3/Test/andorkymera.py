# -*- coding: utf-8 -*-
"""
Created on Mon Nov  6 16:06:04 2023

@author: Lab Di Martino
"""

import nplab
from nplab.instrument.camera.Andor.andor_sdk_rig2 import AndorBase
from nplab.instrument.spectrometer.Kymera import Kymera
from nplab.instrument.spectrometer.seabreeze import OceanOpticsSpectrometer, OceanOpticsControlUI

from nplab.instrument.spectrometer import Spectrometer
import numpy as np

import sys
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap, BoundaryNorm
import matplotlib.colors as colors
import matplotlib.cbook as cbook
from matplotlib.collections import LineCollection



pixel_number = 1024
pixel_width = 26
use_shifts = False
laser_wl = 782.5
white_shutter = None
settings_filepath = None
camera_index = None
parent = None
                     

kymera = Kymera()
kymera.pixel_number = pixel_number
kymera.pixel_width = pixel_width
use_shifts = use_shifts
laser_wl = laser_wl
white_shutter = white_shutter
pixel_number = pixel_number


kymera.SetGrating(grating_num=1)
kymera.SetWavelength(800)

ImageFlip = 0
global wavelength
wavelength = kymera.GetCalibration()[::-1]


""" Initialise the camera """
Andor = AndorBase()
camera_index = None
Andor.start(camera_index)
CurImage = None
background = None
backgrounded = False
keep_shutter_open = False
if settings_filepath is not None:
    load_params_from_file(settings_filepath)
isAborted = False

# %%

available_modes = ['Single', 'Accumulate', 'Kinetic', 'Fast Kinetic']
currentMode = 'Single'
Andor.set_andor_parameter('AcquisitionMode', available_modes.index(currentMode) + 1)
print('Acquisition mode set to Single')

available_modes = ['FVB', 'Multi-track', 'Random track', 'Single track', 'Image']
currentMode = 'Image'
Andor.set_andor_parameter('ReadMode', available_modes.index(currentMode))
print('Read mode set to ' + currentMode)

# %%
Andor.set_andor_parameter('Exposure', float(1))
# %%

imageArray, num_of_images, image_shape = Andor.capture()
imageArray.reverse()


image2D = np.array([ [0] * image_shape[1] for _ in range(int(image_shape[0]))])

for i in np.arange(image_shape[0]):
    for j in np.arange(image_shape[1]):
        image2D[i][j] = imageArray[i * image_shape[1] + j]

        
# %%
%matplotlib qt
real_y = np.arange(image_shape[0])
real_x = wavelength#np.arange(image_shape[1])
real_x.reverse()

real_x = np.arange(image_shape[1])

dx = (real_x[1]-real_x[0])/2.
dy = (real_y[1]-real_y[0])/2.
extent = [real_x[0]-dx, real_x[-1]+dx, real_y[0]-dy, real_y[-1]+dy]

#plt.close('all')
plt.figure(1)
plt.cla()
#clb.remove()
plt.imshow(image2D, extent = extent, aspect='auto', vmin = 280, vmax = 350)
plt.tight_layout()
#clb = plt.colorbar()























