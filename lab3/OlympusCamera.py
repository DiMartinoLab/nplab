# -*- coding: utf-8 -*-
"""
Created on Thu Sep 29 15:43:16 2022

@author: Trung
"""
import nplab
#import nplab
from PyQt5 import QtWidgets, uic
from PyQt5.QtCore import QTimer, QTime, Qt, QRunnable, QThreadPool
import sys
#import os
import numpy as np
import matplotlib.pyplot as plt
#import matplotlib.colors as colors
from nplab.utils.gui import QtCore, QtGui, uic, get_qt_app, show_widget
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg, NavigationToolbar2QT as NavigationToolbar
from nplab.ui.ui_tools import UiTools
import time
from lucam import Lucam
from mpl_toolkits.axes_grid1 import make_axes_locatable


class LiveCapture(QRunnable):
    def __init__(self, Camera_handle, CameraView, cax, Box_ExpTime,\
                 Box_Gain, Box_Brightness, Box_Contrast, Box_Gamma, Box_Saturation):
        
        super().__init__()
        global camera_running
        self.camera = Camera_handle
        self.CameraView = CameraView
        self.cax = cax
        
        self.Box_ExpTime = Box_ExpTime
        self.Box_Gain = Box_Gain
        self.Box_Brightness = Box_Brightness
        self.Box_Contrast = Box_Contrast
        self.Box_Gamma = Box_Gamma
        self.Box_Saturation = Box_Saturation
        image = 0
        image = np.float64(self.camera.TakeSnapshot())
        

        #DFc.SaveImage(image, 'test.tif')
#        plt.figure(1)
#        plt.clf()
#        im = plt.imshow(image/AvFrame - bg, cmap='gray')
#        plt.colorbar(im, orientation='vertical')
#        plt
#        plt.show()

    def run(self):
        global camera_running
        while camera_running == True:
#            print('live')
        
            ExpTime = float(self.Box_ExpTime.value())
            Gain = float(self.Box_Gain.value())
            Brightness = float(self.Box_Brightness.value())
            Contrast = float(self.Box_Contrast.value())
            Gamma = float(self.Box_Gamma.value())
            Saturation = float(self.Box_Saturation.value())
            
            self.camera.set_properties(exposure = ExpTime, gain = Gain,\
                                       contrast = Contrast, brightness = Brightness, \
                                       gamma = Gamma, saturation = Saturation)
            
            image = np.float64(self.camera.TakeSnapshot())
            print(np.amax(image))
            
            start = 10
            stop = 80
            y1 = np.array([-stop, -start]) + len(image)/2
            y11 = np.array([start, stop]) + len(image)/2
            x1 = np.array([len(image[0])/2, len(image[0])/2])
    
            x2 = np.array([-stop, -start]) + len(image[0])/2
            x21 = np.array([start, stop]) + len(image[0])/2
            y2 = np.array([len(image)/2, len(image)/2])
            
            self.CameraView.canvas.ax.cla()
            
    #        self.Temp_Fig.canvas.ax2.cla()
            
            p = self.CameraView.canvas.ax.imshow(image, cmap='binary')
            self.CameraView.canvas.ax.plot(x1, y1, 'r', linewidth = 2)
            self.CameraView.canvas.ax.plot(x1, y11, 'r', linewidth = 2)
            self.CameraView.canvas.ax.plot(x2, y2, 'r', linewidth = 2)
            self.CameraView.canvas.ax.plot(x21, y2, 'r', linewidth = 2)
            
            self.CameraView.canvas.ax.tick_params(axis = "y", direction = "in", left = False, right = False)
            self.CameraView.canvas.ax.tick_params(axis = "x", direction = "in", top = False, bottom = False)
            self.CameraView.canvas.fig.colorbar(p, cax=self.cax)
            self.CameraView.canvas.draw()



        
        
        
        print('Successfully stop')
    
class OlympusCamera(QtWidgets.QWidget, UiTools):
    #inherit QWidget methods
    def __init__(self,  Camera_instance = False, ui_file = 'OlympusCamera.ui',   parent=None):
        
        self.camera = Lucam()
        super(OlympusCamera, self).__init__() 
        uic.loadUi(ui_file, self)
        
        self.live_button.clicked.connect(self.camera_live)
        self.stop_button.clicked.connect(self.camera_stop)
        divider = make_axes_locatable(self.CameraView.canvas.ax)
        self.cax = divider.append_axes("right", size="5%", pad = 0.05) 

        self.camera.set_properties(exposure = 4.5, gain = 0.9, contrast = 68, brightness = 25, gamma = 0.68, saturation = -89)
        self.threadCount = QThreadPool.globalInstance().maxThreadCount()
        self.pool = QThreadPool.globalInstance()
        global camera_running
#        print(str(float(self.Box_Gain.value())))
        print('Initialisation done.')

    def camera_live(self):
        
#        timer = QTimer(self)
#        # Add action to timer
#        timer.timeout.connect(self.temperatureUpdate)
#        # Run the event each 60 seconds
#        timer.start(60*1000)
        global camera_running
        camera_running = True
        camera_live = LiveCapture(self.camera, self.CameraView, self.cax, self.Box_ExpTime, \
                                  self.Box_Gain, self.Box_Brightness, self.Box_Contrast, \
                                  self.Box_Gamma, self.Box_Saturation)
        self.pool.start(camera_live)
        print('bla bla bla stop')
        
    def camera_stop(self):
        global camera_running
        camera_running = False
        image = np.float64(self.camera.TakeSnapshot())
        plt.figure(1)
        plt.clf()
        im = plt.imshow(image, cmap='gray')
        plt.colorbar(im, orientation='vertical')
        
    
        
    def make_window(self):
        app = get_qt_app()
        self.show()
        app.exec_()
        return self

    def plotTemp(self):
        self.CameraView.canvas.ax.cla()
#        self.Temp_Fig.canvas.ax2.cla()
        rate = np.array(self.cMagnetTemp_value[0:len(self.cMagnetTemp_value)-1]) - np.array(self.cMagnetTemp_value[1:len(self.cMagnetTemp_value)])
        
        p = self.CameraView.canvas.ax.plot(np.array(self.Time_count[1:len(self.Time_count)]), rate, 'r')
        self.CameraView.canvas.ax.set_ylabel('Rate (K/m)', color = p[0].get_color())
        self.CameraView.canvas.ax.set_xlabel('Time (m)')
        self.CameraView.canvas.ax.tick_params(axis = "y", direction = "in", left = True, right = False)
        
        
        
        p2 = self.CameraView.canvas.ax2.plot(np.array(self.Time_count), np.array(self.cMagnetTemp_value), 'b')
        self.CameraView.canvas.ax2.set_ylabel('Temperature (K/m)', color = p2[0].get_color())
        self.CameraView.canvas.fig.tight_layout()
        self.CameraView.canvas.ax2.tick_params(axis = "y", direction = "in", left = False, right = True)
        self.CameraView.canvas.ax2.set_ylim([min(np.array(self.cMagnetTemp_value)) - 3, max(np.array(self.cMagnetTemp_value)) + 3])
        self.CameraView.canvas.draw()

app = QtWidgets.QApplication(sys.argv)
main = OlympusCamera()
main.show()
sys.exit(app.exec_())