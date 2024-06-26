
"""
Created on Fri Nov 11 16:35:28 2022

@author: Lab Di Martino
Trung - xtn20
November 2022
Interface and control for the Oxford Mercury iTC system
"""

import nplab
from PyQt5 import QtWidgets, uic
from PyQt5.QtCore import QTimer, QTime, Qt, QRunnable, QThreadPool
import sys
import numpy as np
import os
from nplab.utils.gui import QtCore, QtGui, uic, get_qt_app, show_widget
from nplab.ui.ui_tools import UiTools
from nplab.instrument.mercuryUSB.mercuryUSB import temperatureController as MiC_package
import time
from datetime import datetime


#import curve to be used when using reflective metal as reference material rather than white scatterer
ratio_curve_points = np.genfromtxt(fname='ratio_curve.txt')

class sub_read_current_temps(QRunnable):
    def __init__(self, MiTC_handle, cSampleTempDisplay, cMagnetTempDisplay, cHTSTempDisplay):
        super().__init__()
        self.MiTC_handle = MiTC_handle
        # 3 display fields on the UI for sample, magnet and HTS
        self.cSampleTempDisplay = cSampleTempDisplay
        self.cMagnetTempDisplay = cMagnetTempDisplay
        self.cHTSTempDisplay = cHTSTempDisplay
        self.cTime_value = 0
        self.cSampleTemp_value = 0
        self.cMagnetTemp_value = 0
        
    
    def run(self):
        # check and block the communication with MiTC while reading
        global MiTC_communication
        global reading_done
        # check the communication is free, if not, wait till it is freed
        while MiTC_communication > 0:
            time.sleep(0.05)
        # block the communication
        MiTC_communication = 1
        # read temperature here
        self.cTime_value = str(self.MiTC_handle.getTime())
        self.cSampleTemp_value = str(self.MiTC_handle.getTemp('Sample'))
        self.cSampleTempDisplay.display(float(self.cSampleTemp_value))
        self.cMagnetTemp_value = str(self.MiTC_handle.getTemp('Magnet'))
        self.cMagnetTempDisplay.display(float(self.cMagnetTemp_value))
        self.cHTSTemp_value = str(self.MiTC_handle.getTemp('HTS'))
        self.cHTSTempDisplay.display(float(self.cHTSTemp_value))
        # free the communication
        MiTC_communication = 0
        reading_done = 1

class MercuryController(QtWidgets.QWidget, UiTools):
    #inherit QWidget methods
    def __init__(self, MiTC_handle = False,  ui_file = 'MercuryController.ui',   parent=None):
        
#        # set objects for spectrometer and stage
#        self.MiTC_handle = MiTC_handle
        print('Mercury controller widget has been requested')
        self.MiTC_handle = MiC_package('COM8')
        
        # define devices of MiTC and MiPS.
        # For more information, open "mercurytest.py" in
        # nplab/instrument/mercuryUSB
        self.MiTC_handle.devices = {'Magnet' : 'DEV:DB1.T1:TEMP', 'HTS' : 'DEV:DB2.T1:TEMP', 'Sample' : 'DEV:MB1.T1:TEMP', 'Heater' : 'DEV:MB0.H1:HTR'}
        # make this window available for the main program
        super(MercuryController, self).__init__() 
        uic.loadUi(ui_file, self)
        """ """
        now = datetime.now()
        current_date = now.strftime("%y%m%d")
        
        
        self.folder2save = os.path.join('C:\\Users\\Lab Di Martino\\Documents\\data\\xtn20')
        
        folders = []
        for root, dirs, files in os.walk(self.folder2save):
        
            for folder in dirs:
                if folder.startswith(current_date):
                    folders.append(folder)
        if not os.path.exists(self.folder2save + '\\' + str(folders[0])):
            os.makedirs(self.folder2save + '\\' + str(folders[0]))
            
        path2save = os.path.join(self.folder2save + '\\' + str(folders[0]) + '\\Cooldown')
        if not os.path.exists(path2save):
            os.makedirs(path2save)
            
        """ """
        self.folder = path2save
        """ """
         # Save as increment file name
        i = 0
        
        fullfile = path2save + '\\Templog_TSMH'
        if os.path.exists(fullfile + '.txt'):
            while os.path.exists(fullfile + '_' + str(i) + '.txt'):
                print('file found')
                i = i + 1
            print(i)
            outputfile = fullfile + '_' + str(i) + '.txt'
        else:
            outputfile = fullfile + '.txt'
            
       
        """ """
        self.file_name = outputfile
        self.cTime_value = []
        self.Time_count = []
        self.cSampleTemp_value = []
        self.cMagnetTemp_value = []
        self.cHTSTemp_value = []
        self.count = 0
        self.Start_button.clicked.connect(self.startTimer)
        self.threadCount = QThreadPool.globalInstance().maxThreadCount()
        self.pool = QThreadPool.globalInstance()
        
        self.test_button.clicked.connect(self.test_function)

    def test_function(self):
        self.MiTC_handle.setSampleTemp('Sample', 10)
        
    def startTimer(self):
        print('Start logging')
        self.temperatureUpdate()
        self.plotTemp()
        # Create a timer object
        timer = QTimer(self)
        # Add action to timer
        timer.timeout.connect(self.temperatureUpdate)
        # Run the event each 60 seconds
        timer.start(60*1000)

        
    def temperatureUpdate(self):
        read_temp = sub_read_current_temps(self.MiTC_handle, self.cSampleTemp, self.cMagnetTemp, self.cHTSTemp)
        
        global MiTC_communication
        global reading_done
        reading_done = 0
        MiTC_communication = 0
        self.pool.start(read_temp)
        
        while reading_done < 1:
            test = 0
        reading_done = 0
        self.cTime_value.append((read_temp.cTime_value))
        self.cSampleTemp_value.append(float(read_temp.cSampleTemp_value))
        self.cMagnetTemp_value.append(float(read_temp.cMagnetTemp_value))
        self.cHTSTemp_value.append(float(read_temp.cHTSTemp_value))
        self.Time_count.append(self.count)
        self.count = self.count + 1
        self.plotTemp()
        savestr = str(read_temp.cTime_value) + '\t' + str(read_temp.cSampleTemp_value) \
            + '\t' + str(read_temp.cMagnetTemp_value)  + '\t' + str(read_temp.cHTSTemp_value) + '\n'
        with open(self.file_name, 'a') as the_file:
            the_file.write(savestr)
        
    def plotTemp(self):
        self.Temp_Fig.canvas.ax.cla()
#        self.Temp_Fig.canvas.ax2.cla()
        rate = np.array(self.cMagnetTemp_value[0:len(self.cMagnetTemp_value)-1]) - np.array(self.cMagnetTemp_value[1:len(self.cMagnetTemp_value)])
        
        p = self.Temp_Fig.canvas.ax.plot(np.array(self.Time_count[1:len(self.Time_count)]), rate, 'r')
        self.Temp_Fig.canvas.ax.set_ylabel('Rate (K/m)', color = p[0].get_color())
        self.Temp_Fig.canvas.ax.set_xlabel('Time (m)')
        self.Temp_Fig.canvas.ax.tick_params(axis = "y", direction = "in", left = True, right = False)
        
        minY = min(float(self.Rate1.value()), float(self.Rate2.value()))
        maxY = max(float(self.Rate1.value()), float(self.Rate2.value()))
        if minY<maxY:
            self.Temp_Fig.canvas.ax.set_ylim([minY, maxY])
        
        minX = min(float(self.Time1.value()), float(self.Time2.value()))
        maxX = max(float(self.Time1.value()), float(self.Time2.value()))
        if minX<maxX:
            self.Temp_Fig.canvas.ax.set_xlim([minX, maxX])
        
        p2 = self.Temp_Fig.canvas.ax2.plot(np.array(self.Time_count), np.array(self.cMagnetTemp_value), 'b')
        self.Temp_Fig.canvas.ax2.plot(np.array(self.Time_count), np.array(self.cSampleTemp_value), 'b+')
        self.Temp_Fig.canvas.ax2.plot(np.array(self.Time_count), np.array(self.cHTSTemp_value), 'bo')
        self.Temp_Fig.canvas.ax2.set_ylabel('Temperature (K/m)', color = p2[0].get_color())
        self.Temp_Fig.canvas.fig.tight_layout()
        self.Temp_Fig.canvas.ax2.tick_params(axis = "y", direction = "in", left = False, right = True)
        if minX<maxX:
            self.Temp_Fig.canvas.ax2.set_xlim([minX, maxX])
            
        minY = min(float(self.Temp1.value()), float(self.Temp2.value()))
        maxY = max(float(self.Temp1.value()), float(self.Temp2.value()))
        if minY<maxY:
            self.Temp_Fig.canvas.ax2.set_ylim([minY, maxY])
        else:
            self.Temp_Fig.canvas.ax2.set_ylim([min(np.array(self.cMagnetTemp_value)) - 3, max(np.array(self.cMagnetTemp_value)) + 3])
        self.Temp_Fig.canvas.draw()
        
        
    
#    def make_window(self):
#        app = get_qt_app()
#        self.show()
#        app.exec_()
#        return self

app = QtWidgets.QApplication(sys.argv)
main = MercuryController()
main.show()
sys.exit(app.exec_())











