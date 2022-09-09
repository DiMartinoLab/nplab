# -*- coding: utf-8 -*-
"""
Created on Wed Jul 20 11:13:29 2022

@author: Lab Di Martino
"""

"""lab3_interface_dawn_troubleshooting"""

# -*- coding: utf-8 -*-
"""
Created on Jan 15 10:23:36 2019

@author: Hera
"""

import nplab

from nplab.instrument.spectrometer.shamdor import Shamdor
from nplab.instrument.spectrometer.seabreeze import OceanOpticsSpectrometer
from nplab.instrument.stage.SMC100 import SMC100
import nplab.utils.gui 
import nplab.datafile as datafile
from nplab.instrument.spectrometer import Spectrometer
from nplab.experiment import Experiment, ExperimentStopped
from nplab.utils.notified_property import DumbNotifiedProperty, NotifiedProperty
from nplab.utils.gui import QtCore, QtGui, uic, get_qt_app
from nplab.ui.ui_tools import UiTools
import os
import time
import threading
import numpy as np
import pyqtgraph as pg
from nplab.ui.ui_tools import QtWidgets



#myOceanOptics = OceanOpticsSpectrometer(0)

class Lab3_experiment(Experiment, QtWidgets.QWidget, UiTools):
    # To use auto_connect_by_name name all widgets using _WidgetType, e.g. Vhigh_DoubleSpinBox
    # Then define DumbNotifiedProperty with the same name without _WidgetType, e.g. Vhigh

    smu_wait = DumbNotifiedProperty(0.0)
    
    centre_row = DumbNotifiedProperty(100)
    num_rows = DumbNotifiedProperty(15)

    description = DumbNotifiedProperty("description")
    log_to_console = True
    live_darkfield_spectrum_signal = QtCore.Signal(np.ndarray)

    
    def __init__ (self, activeDatafile):
        super(Lab3_experiment, self).__init__()
#        uic.loadUi('lab3_interface.ui', self)
        uic.loadUi('lab3_interface_sunny.ui', self)
        
###comment out software you are not going to use

        self.initialise_SMC100() #actuators for xy stage
        self.initialise_OOSpectrometer() #for DF (white light) and PL (444nm laser)

####end        
        self.setup_plot_widgets()
        
        self.activeDatafile = activeDatafile

        self.auto_connect_by_name(self)
        
        self.live_darkfield_spectrum_signal.connect(self.update_darkfield_spectrum_plot)
        self.openstage.clicked.connect(self.open_SMC100_ui)

    def run(self, *args):
        # *args collects extra unnecessary arguments from qt

        self.times_data = []

        stepwiseCounter = 1
        rampCounter = 1
        runningRampInterval = False
        try:
            activeDatagroup = self.create_data_group('scan_%d')
            activeDatagroup.attrs.create('description', str(self.description))
            if not self.mode_smuOnly.isChecked():
                
                if (self.mode_DarkfieldOnly.isChecked() or self.mode_RamanAndDarkfield.isChecked() ):
                    self.darkfieldSpectrumImagePlotData = []
                    self.darkfieldWavelengths = self.OOspectrometer.read_wavelengths()
                    activeDatagroup.attrs.create('DarkfieldWavelengths', self.darkfieldWavelengths)
                    activeDatagroup.attrs.create('DarkfieldBackground', self.OOspectrometer.background)
                    activeDatagroup.attrs.create('DarkfieldReference', self.OOspectrometer.reference)
                    activeDatagroup.attrs.create('DarkfieldIntegrationTime', self.OOspectrometer.get_integration_time())
                    activeDatagroup.attrs.create('DarkfieldBackground_Reference', self.OOspectrometer.background_ref)
                    activeDatagroup.attrs.create('DarkfieldBackground_ReferenceIntegrationTime', self.OOspectrometer.background_int_ref)                    
                    #TEST BY THOMAS
                    #adding attributes for DF in-situ measurements
                    activeDatagroup.attrs.create('DarkfieldReferenceIntegrationTime', self.OOspectrometer.reference_int)
                    activeDatagroup.attrs.create('DarkfieldBackgroundIntegrationTime', self.OOspectrometer.background_int)
                if (self.mode_PLOnly.isChecked()):
                    self.darkfieldSpectrumImagePlotData = []
                    self.darkfieldWavelengths = self.OOspectrometer.read_wavelengths()
                    activeDatagroup.attrs.create('PLWavelengths', self.darkfieldWavelengths)
                    activeDatagroup.attrs.create('PLBackground', self.OOspectrometer.background)
                    activeDatagroup.attrs.create('PLIntegrationTime', self.OOspectrometer.get_integration_time())
                    
        except ExperimentStopped:
            print("----- Measurement stopped -----")
        finally:
            pass

           
    def initialise_SMC100(self):
        self.SMC100=SMC100('COM1', (1,2,3))   
        
    def initialise_OOSpectrometer(self):
        self.OOspectrometer = OceanOpticsSpectrometer(0)

    
    def setup_plot_widgets(self):
        print('In development')

        self.darkfieldSpectrum_plot = pg.PlotWidget()
        self.darkfieldSpectrum_vs_time_plot = pg.PlotWidget()
        self.replace_widget(self.plotGrid, self.plot6, self.darkfieldSpectrum_vs_time_plot)
#        # have to use pg.ImageItem() for an image plot. ImageItem is not a widget so it can't replace a PlotWidget directly,
#        # but it can be added as an ImageItem to a PlotWidget
        self.darkfieldImagePlotItem = pg.ImageItem()
        self.darkfieldSpectrum_vs_time_plot.addItem(self.darkfieldImagePlotItem)
        

    def update_darkfield_spectrum_plot(self, spectrum):
        self.darkfieldSpectrum_plot.plot(self.darkfieldWavelengths, spectrum, clear = True, pen = 'r')  # plot current darkfield spectrum in real time
        self.darkfieldSpectrumImagePlotData.append(np.nan_to_num(spectrum)[80:750])     # plot darkfield-time in range 400-900nm as image plot. ImageItem doesn't handle nan values, converting those to zero   
        self.darkfieldImagePlotItem.setImage(np.asarray(self.darkfieldSpectrumImagePlotData))
        
        
    def open_SMC100_ui(self):
        self.SMC100.show_gui()
        
   
    def acquire_darkfield_spectrum(self):
        self.darkfieldSpectrum = np.asarray(self.OOspectrometer.read_processed_spectrum())
        
    def open_OO_spectrometer(self):
#        self.OOspectrometer = OceanOpticsSpectrometer(0)
        self.gui_OOspectrometer= self.OOspectrometer.get_qt_ui()
        self.gui_OOspectrometer.show()
    
    def shutdown(self):
        self.activeDatafile.close()
        self.myShamdor.cooler = False
        self.myShamdor.close()
        self.OOspectrometer.shutdown_seabreeze()
        print('----Experiment ended----')
        #TODO - closes any UI windows
    
    def get_qt_ui(self):
        return self
    

if __name__ == '__main__':
    activeDatafile = nplab.current_datafile(working_directory = os.path.abspath(os.path.join(os.getcwd(),"../../..")))
    #working directory should now be C:\\Users\\<name>\\Documents    
    gui_activeDatafile = activeDatafile.get_qt_ui()
    gui_activeDatafile.show()
    
    experiment = Lab3_experiment(activeDatafile)
    experiment.show_gui()
    print('Done')
