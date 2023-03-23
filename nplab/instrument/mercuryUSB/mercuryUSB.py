"""
Created on Fri Nov 11 16:35:28 2022

@author: Lab Di Martino
Trung - xtn20
November 2022 - Control and calibration routine for the Oxford Mercury iTC and iPS-M system

"""

import serial
import numpy as np
import time
import socket
import sys

class temperatureController(object):
    def __init__(self, port, baudrate = 115200):
        self.ser = serial.Serial(port, baudrate = baudrate, stopbits = serial.STOPBITS_ONE, timeout = 1)

    def defineDevices(self, devices_input):
        #store the device addresses in a dictionary
        self.devices = devices_input

    def getVersion(self):
        string = self.readValue("*IDN?", readPrefix = "")
        return string

    def getDevices(self):
        devices = self.readValue("SYS:CAT")
        return devices

    def writeValue(self, value):
        self.ser.write(bytes(value + "\n", 'UTF-8'))
        time.sleep(3)

    def readValue(self, value, readPrefix = "READ:"):
        self.writeValue(readPrefix + value)
        string = self.ser.readline().rstrip()
        self.ser.flush()
        return string

    def setValue(self, value):
        self.writeValue("SET:" + value)
        string = self.ser.readline().rstrip()

    def close(self):
        self.ser.close()

    def getTime(self):
        ctime = str(self.readValue('SYS:TIME'))
        ctime = ctime[len(ctime)-9:len(ctime)-1]        
        return ctime

    def getSignal(self, device, signal):
        """Get a signal from a device
        unitPrefixes are taken into account and the value is returned as a float.

        - device: device key for the devices dictionary
        - signal: string corresponding to a valid signal, i.e. TEMP, VOLT, CURR, RES, etc.
        """
        ans = str(self.readValue(self.devices[device] + ":SIG:" + signal)).split(":")[-1]
        
        siPrefixes = {"M": 1e6, "k" : 1e3, "m" : 1e-3, "\xb5" : 1e-6, "n" : 1e-9, "p" : 1e-12}
        
        if ans[-2].isdigit():
            return float(ans[:-2])
        else:
            try:
                return float(ans[:-2])*siPrefixes[ans[-2]]
            except:
                print( "Ans: ", ans)
                raise ValueError 
        
    def getTemp(self, device):
        ans = str(self.readValue(self.devices[device] + ":SIG:TEMP")).split(":")[-1]
        return float(ans[0:len(ans)-2])

    def getField(self, device):
        ans = str(self.readValue(self.devices[device] + ":SIG:FLD")).split(":")[-1]
        return float(ans[0:len(ans)-2])
    
    def getFieldRate(self, device, value):
        ans = str(self.readValue(self.devices[device] + ":SIG:RFST:" + str(value))).split(":")[-1]
        print(ans[0:len(ans)-2])
    
    def getMiPSHeater(self, device):
        ans = str(self.readValue(self.devices[device] + ":SIG:SWHT")).split(":")[-1]
        return float(ans[0:len(ans)-2])
    
    def setMiPSHeater(self, device, state):
        ans = str(self.setValue(self.devices[device] + ":SIG:SWHT:" + state)).split(":")[-1]
        return float(ans[0:len(ans)-2])
    
    def setField(self, device, value):
        self.setValue(self.devices[device] + ":SIG:FSET:" + str(value))

    def setFieldRate(self, device, value):
        self.setValue(self.devices[device] + ":SIG:RFST:" + str(value))

    def setAction(self, device, action):
        """action = [HOLD | RTOS | RTOZ | CLMP]
        RTOS= Ramp to set         RTOZ= Ramp to zer       CLMP= clamp output if CURR=0."""
        self.setValue(self.devices[device] + ":ACTN:" + action)
        print(self.devices[device] + ":ACTN:" + action)
     
    
    def getSensorInformation(self, device, includeTemperature = False):
        """Get Voltage, Current, Resistance and optionally Temperature of a device"""

        v = self.getSignal(device, "VOLT")
        c = self.getSignal(device, "CURR")
        r = self.getSignal(device, "RES")

        if includeTemperature:
            t = self.getSignal(device, "TEMP")
            return v,c,r,t
        else:
            return v,c,r
    
    def calibrate(self, sensorNameDB6, sensorNameMB):
        """Calibration Routine. Calibrated sensor connected to port DB7
        Two sensors to be calibrated connected to ports DB6 and MB
        
        The temperature, resistance, current and voltage of the calibrated sensor is recorded first, then current and voltage of the second sensor are recorded, then current and voltage of the last sensor are recorded.

        Make sure that all boards are configured correctly, i.e. for NTC sensors in constant voltage mode, with voltage set to 7 mV. Configuration has to be done using the devices Touchscreen UI.

        Note that this routine only records the three sensor values. In our setup the temperature sweep is achieved using a second temperature controller.
        """

        #in the command below, the unit of magnitude is mV and should not be supplied.
        #everything is working except for the calibration file
        #This currently has to be set via the Touchscreen Interface of the device.
        #self.setValue(self.devices["db7"] + ":TYPE:NTC:EXCT:TYPE:UNIP:MAG:7:CALB:X96620.dat")

        calList = []
        sensor_db6List = []
        sensor_mb1List = []
        
        while True:
            try:
                cal = self.getSensorInformation("db7", includeTemperature = True)
                sensor_db6 = self.getSensorInformation("db6")
                sensor_mb1 = self.getSensorInformation("mb1")

                calList.append(cal)
                sensor_db6List.append(sensor_db6)
                sensor_mb1List.append(sensor_mb1)

                exportArray = np.hstack((np.array(calList), np.array(sensor_db6List), np.array(sensor_mb1List)))

                print( "Export Array, ", exportArray)
                print( "Rc: ", cal[2], "T: ", cal[3], "R1: ", sensor_db6[2], "R2: ", sensor_mb1[2])
                
                headerString = """######################################

  Calibration Log

  Calibrated sensor: X_____.dat 
  (Columns 1 to 4: Voltage, Current, Resitance, Temperature)

  First sensor to calibrate:""" + sensorNameDB6 + """
  (Columns 5 to 7: Voltage, Current, Resistance)

  Second sensor to calibrate:""" + sensorNameMB + """
  (Columns 8 to 10: Voltage, Current, Resistance)

######################################"""
                
                np.savetxt("calibration.txt", exportArray, fmt="%.6e", header = headerString, comments = "#")

                np.savetxt(sensorNameDB6 + ".dat", exportArray[:,[3,6]], fmt="%.6e", header = "Temperature (K)\t Resistance (Ohm)\nExcitation: Constant Voltage, 7mV", comments = "#")
                np.savetxt(sensorNameMB + ".dat", exportArray[:,[3,9]], fmt="%.6e", header = "Temperature (K)\t Resistance (Ohm)\nExcitation: Constant Voltage, 7mV", comments = "#")
                time.sleep(1)        
                                
            except KeyboardInterrupt:
                print( "Keyboard Interrupt caught. Finishing Calibration.")
                print( "Minimum Temperature achieved: ", np.min(exportArray[:,3]))
                print( "Maximum Temperature achieved: ", np.max(exportArray[:,3]))
                break

    def autoPollTemperatures(self):
        timeList = []
        t1List = []
        t2List = []
        t3List = []

        startTime = time.time()
        while True:
            try:
                timeNow = time.time() - startTime
                t1 = self.getSignal("mb1", "TEMP")
                t2 = self.getSignal("db6", "TEMP")
                t3 = self.getSignal("db7", "TEMP")

                print( "MB1 {:.3f} K    DB6 {:.3f} K       DB7 {:.3f} K".format(t1, t2, t3))

                timeList.append(timeNow)
                t1List.append(t1)
                t2List.append(t2)
                t3List.append(t3)

                #np.savetxt("tempLog_" + str(startTime) + ".txt", np.transpose(np.vstack((np.array(timeList), np.array(t1List), np.array(t2List), np.array(t3List)))), header = "Time\tT1\tT2\tT3", comments="#", fmt="%.6e")

                np.savetxt("tempLog_" + str(startTime) + ".txt", np.vstack((np.array(timeList), np.array(t1List), np.array(t2List), np.array(t3List))).T, header = "Time\tT1\tT2\tT3", comments="#", fmt="%.6e")

                time.sleep(1)
            except KeyboardInterrupt:
                break


            