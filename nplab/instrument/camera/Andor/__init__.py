import nplab.utils.gui
from nplab.utils.gui import QtWidgets, QtCore, uic
from nplab.instrument.camera import Camera
from nplab.utils.thread_utils import background_action, locked_action

import os
import platform
import sys
import time
from ctypes import *
import numpy as np
from PIL import Image


class AndorCapabilities(Structure):
    _fields_ = [("ulSize", c_ulong),
                ("ulAcqModes", c_ulong),
                ("ulReadModes", c_ulong),
                ("ulTriggerModes", c_ulong),
                ("ulCameraType", c_ulong),
                ("ulPixelMode", c_ulong),
                ("ulSetFunctions", c_ulong),
                ("ulGetFunctions", c_ulong),
                ("ulFeatures", c_ulong),
                ("ulPCICard", c_ulong),
                ("ulEMGainCapability", c_ulong),
                ("ulFTReadModes", c_ulong)]


class AndorWarning(Warning):
    def __init__(self, code, msg, reply):
        super(AndorWarning, self).__init__()
        self.error_code = code
        self.error_name = ERROR_CODE[code]

        self.msg = msg
        self.reply = reply

    def __str__(self):
        return self.error_name + '\n Error sent: ' + self.msg + '\n Error reply: ' + self.reply


class AndorBase:
    """
    The self.parameters dictionary contains all the information necessary to deal with the camera parameters. Each
    entry in the dictionary corresponds to a specific parameter and allows you to specify the Get and/or Set command
    name and datatype (from the .dll).

    Most parameters are straightforward, since the Andor dll either has inputs (for setting parameters) or outputs
    (for getting parameters). So you can just intuitively call Andor.GetParameter(name) or Andor.SetParameter(name, value)
    with name and value provided by the user.
    Some parameters, like VSSpeed, HSSpeed..., require inputs to get outputs, so the user must say, e.g.,
        Andor.GetParameter('VSSpeed', 0)
    Which does not return the current VSSpeed, but the VSSpeed (in us) of the setting 0.
    """

    def __init__(self):
        if platform.system() == 'Windows':
            if platform.architecture()[0] == '32bit':
                self.dll = windll(os.path.dirname(__file__) + "\\atmcd32d")
            elif platform.architecture()[0] == '64bit':
                self.dll = CDLL(os.path.dirname(__file__) + "\\atmcd64d")
            else:
                raise Exception("Cannot detect Windows architecture")
        elif platform.system() == "Linux":
            dllname = "usr/local/lib/libandor.so"
            self.dll = cdll.LoadLibrary(dllname)
        else:
            raise Exception("Cannot detect operating system for Andor")

        self.parameters = dict(
            SoftwareWaitBetweenCaptures=(),
            DetectorShape=dict(Get=dict(cmdName='GetDetector', Outputs=(c_int, c_int)), value=None),
            SerialNumber=dict(Get=dict(cmdName='GetCameraSerialNumber', Outputs=(c_int,)), value=None),
            HeadModel=dict(Get=dict(cmdName='GetHeadModel', Outputs=(c_char,) * 20), value=None),
            Capabilities=dict(Get=dict(cmdName='GetCapabilities', Outputs=(
                              AndorCapabilities(sizeof(c_ulong) * 12, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),)), value=None),
            AcquisitionMode=dict(Set=dict(cmdName='SetAcquisitionMode', Inputs=(c_int,)), value=None),
            TriggerMode=dict(Set=dict(cmdName='SetTriggerMode', Inputs=(c_int,)), value=None),
            ReadMode=dict(Set=dict(cmdName='SetReadMode', Inputs=(c_int,)), value=None),
            CropMode=dict(Set=dict(cmdName='SetCropMode', Inputs=(c_int,) * 3), value=None),
            IsolatedCropMode=dict(Set=dict(cmdName='SetIsolatedCropMode', Inputs=(c_int,) * 5), value=(0,)),
            AcquisitionTimings=dict(Get=dict(cmdName='GetAcquisitionTimings', Outputs=(c_float, c_float, c_float)),
                                    value=None),
            AccumCycleTime=dict(Set=dict(cmdName='SetAccumulationCycleTime', Inputs=(c_float,)),
                                Finally='AcquisitionTimings'),
            KinCycleTime=dict(Set=dict(cmdName='SetKineticCycleTime', Inputs=(c_float,)),
                              Finally='AcquisitionTimings'),
            Exposure=dict(Set=dict(cmdName='SetExposureTime', Inputs=(c_float,)), Finally='AcquisitionTimings'),
            Image=dict(Set=dict(cmdName='SetImage', Inputs=(c_int,) * 6), value=None),
            NAccum=dict(Set=dict(cmdName='SetNumberAccumulations', Inputs=(c_int,)), value=1),
            NKin=dict(Set=dict(cmdName='SetNumberKinetics', Inputs=(c_int,)), value=1),
            FastKinetics=dict(Set=dict(cmdName='SetFastKineticsEx', Inputs=(c_int, c_int, c_float,) + (c_int,) * 4),
                              value=None),
            EMGain=dict(Set=dict(cmdName='SetEMCCDGain', Inputs=(c_int,)),
                        Get=dict(cmdName='GetEMCCDGain', Outputs=(c_int,)), value=None),
            EMAdvancedGain=dict(Set=dict(cmdName='SetEMAdvanced', Inputs=(c_int,)), value=None),
            EMMode=dict(Set=dict(cmdName='SetEMCCDGainMode', Inputs=(c_int,)), value=None),
            EMGainRange=dict(Set=dict(cmdName='GetEMCCDGainRange', Outputs=(c_int,) * 2), value=None),
            Shutter=dict(Set=dict(cmdName='SetShutter', Inputs=(c_int,) * 4), value=None),
            CoolerMode=dict(Set=dict(cmdName='SetCoolerMode', Inputs=(c_int,)), value=None),
            FanMode=dict(Set=dict(cmdName='SetFanMode', Inputs=(c_int,)), value=None),
            ImageFlip=dict(Set=dict(cmdName='SetImageFlip', Inputs=(c_int,) * 2), value=None),
            ImageRotate=dict(Set=dict(cmdName='SetImageRotate', Inputs=(c_int,)), value=None),
            CurrentTemperature=dict(Get=dict(cmdName='GetTemperature', Outputs=(c_int,)), value=None),
            SetTemperature=dict(Set=dict(cmdName='SetTemperature', Inputs=(c_int,)), value=None),
            OutAmp=dict(Set=dict(cmdName='SetOutputAmplifier', Inputs=(c_int,))),
            FrameTransferMode=dict(Set=dict(cmdName='SetFrameTransferMode', Inputs=(c_int,)), value=None),
            SingleTrack=dict(Set=dict(cmdName='SetSingleTrack', Inputs=(c_int,) * 2), value=None),
            MultiTrack=dict(Set=dict(cmdName='SetMultiTrack', Inputs=(c_int,) * 3, Outputs=(c_int,) * 2)),
            FVBHBin=dict(Set=dict(cmdName='SetFVBHBin', Inputs=(c_int,)), value=None),
            Spool=dict(Set=dict(cmdName='SetSpool', Inputs=(c_int, c_int, c_char, c_int)), value=None),
            NumVSSpeed=dict(Get=dict(cmdName='GetNumberVSSpeeds', Outputs=(c_int,)), value=None),
            NumHSSpeed=dict(Get=dict(cmdName='GetNumberHSSpeeds', Inputs=(c_int, c_int,), Outputs=(c_int,)),
                            value=None),
            VSSpeed=dict(Set=dict(cmdName='SetVSSpeed', Inputs=(c_int,)),
                         Get=dict(cmdName='GetVSSpeed', Inputs=(c_int,), Outputs=(c_float,)), GetAfterSet=True),
            HSSpeed=dict(Set=dict(cmdName='SetHSSpeed', Inputs=(c_int, c_int)),
                         Get=dict(cmdName='GetHSSpeed', Inputs=(c_int,) * 3, Outputs=(c_float,))),
            NumPreAmp=dict(Get=dict(cmdName='GetNumberPreAmpGains', Outputs=(c_int,))),
            PreAmpGain=dict(Set=dict(cmdName='SetPreAmpGain', Inputs=(c_int,)),
                            Get=dict(cmdName='GetPreAmpGain', Inputs=(c_int,), Outputs=(c_float,)), GetAfterSet=True),
            NumADChannels=dict(Get=dict(cmdName='GetNumberADChannels', Outputs=(c_int,))),
            ADChannel=dict(Set=dict(cmdName='SetADChannel', Inputs=(c_int,))),
            BitDepth=dict(Get=dict(cmdName='GetBitDepth', Inputs=(c_int,), Outputs=(c_int,)))
        )

        self.Initialize()
        self.capabilities = {}

    def __del__(self):
        """
        If the camera is not a Newton, we start a thread that determines whether the camera is too cold for shutdown.
        This can be modified to include other cameras that are safe to shutdown when cold
        :return:
        """
        if self.parameters['Capabilities']['value']['CameraType'] != 8:
            waitThread = WaitThread(self)
            waitThread.start()
            waitThread.wait()
        self._dllWrapper('ShutDown')
        libHandle = self.dll._handle
        if platform.system() == 'Windows':
            if platform.architecture()[0] == '32bit':
                windll.kernel32.FreeLibrary(libHandle)
            elif platform.architecture()[0] == '64bit':
                # Following http://stackoverflow.com/questions/19547084/can-i-explicitly-close-a-ctypes-cdll
                from ctypes import wintypes
                kernel32 = WinDLL('kernel32')
                kernel32.FreeLibrary.argtypes = [wintypes.HMODULE]
                kernel32.FreeLibrary(libHandle)
            else:
                raise Exception("Cannot detect Windows architecture")
        elif platform.system() == "Linux":
            cdll.LoadLibrary('libdl.so').dlclose(libHandle)
        del self.dll

    '''Base functions'''

    @locked_action
    def _dllWrapper(self, funcname, inputs=(), outputs=(), reverse=False):
        """Handler for all the .dll calls of the Andor

        Parameters
        ----------
        funcname    Name of the dll function to be called
        inputs      Inputs to be handed in to the dll function
        outputs     Outputs to be expected from the dll
        reverse     Whether to have the inputs first or the outputs first when calling the dll

        Returns
        -------

        """
        dll_input = ()
        if reverse:
            for output in outputs:
                dll_input += (byref(output),)
            for inpt in inputs:
                dll_input += (inpt['type'](inpt['value']),)
        else:
            for inpt in inputs:
                dll_input += (inpt['type'](inpt['value']),)
            for output in outputs:
                dll_input += (byref(output),)

        error = getattr(self.dll, funcname)(*dll_input)
        self._errorHandler(error, funcname, *(inputs + outputs))

        returnVals = ()
        for output in outputs:
            if hasattr(output, 'value'):
                returnVals += (output.value,)
            if isinstance(output, AndorCapabilities):
                dicc = {}
                for key, value in output._fields_:
                    dicc[key[2:]] = getattr(output, key)
                returnVals += (dicc,)
        if len(returnVals) == 1:
            return returnVals[0]
        else:
            return returnVals

    def _errorHandler(self, error, funcname='', *args):
        if '_logger' in self.__dict__:
            self._logger.debug("[%s]: %s %s" % (funcname, ERROR_CODE[error], str(args)))
        elif 'verbosity' in self.__dict__:
            if self.verbosity:
                print "[%s]: %s" % (funcname, ERROR_CODE[error])
        if funcname == 'GetTemperature':
            return
        if error != 20002:
            raise AndorWarning(error, funcname, ERROR_CODE[error])

    @background_action
    def _constantlyUpdateTemperature(self):
        self.aborted = False
        while not self.aborted:
            print self.GetParameter('CurrentTemperature')
            time.sleep(10)

    def SetParameter(self, param_loc, *inputs):
        """Parameter setter

        Using the information contained in the self.parameters dictionary, send a general parameter set command to the
        Andor. The command name, and number of inputs and their types are stored in the self.parameters

        Parameters
        ----------
        param_loc   dictionary key of self.parameters
        inputs      inputs required to set the particular parameter. Must be at least one

        Returns
        -------

        """
        func = self.parameters[param_loc]['Set']

        form_in = ()
        for ii in range(len(inputs)):
            form_in += ({'value': inputs[ii], 'type': func['Inputs'][ii]},)

        self._dllWrapper(func['cmdName'], inputs=form_in)

        if len(inputs) == 1:
            self.parameters[param_loc]['value'] = inputs[0]
        else:
            self.parameters[param_loc]['value'] = inputs

        if 'Finally' in self.parameters[param_loc]:
            self.GetParameter(self.parameters[param_loc]['Finally'])
        if 'GetAfterSet' in self.parameters[param_loc]:
            self.GetParameter(param_loc, *inputs)

    def GetParameter(self, param_loc, *inputs):
        """Parameter getter

        Using the information contained in the self.parameters dictionary, send a general parameter get command to the
        Andor. The command name, and number of inputs and their types are stored in the self.parameters

        Parameters
        ----------
        param_loc   dictionary key of self.parameters
        inputs      optional inputs for getting the specific parameter

        Returns
        -------

        """
        if 'Get' in self.parameters[param_loc].keys():
            func = self.parameters[param_loc]['Get']

            form_out = ()
            if param_loc == 'Capabilities':
                form_out += (func['Outputs'][0],)
            else:
                for output in func['Outputs']:
                    form_out += (output(),)
            form_in = ()
            for ii in range(len(inputs)):
                form_in += ({'value': inputs[ii], 'type': func['Inputs'][ii]},)

            vals = self._dllWrapper(func['cmdName'], inputs=form_in, outputs=form_out)
            # if len(vals) == 1:
            #     vals = vals[0]
            self.parameters[param_loc]['value'] = vals
            return vals
        else:
            return self.parameters[param_loc]['value']

    def GetAllParameters(self):
        '''Gets all the parameters that can be gotten

        The parameters that can be gotten are those that have a get capability in the .dll
        Parameter getters that require inputs (i.e. HSSpeed, VSSpeed, PreAmpGain, BitDepth and NumHSSpeed), have to be
        handled separately, and in particular HSSpeed, VSSpeed and PreAmpGain cannot be retrieved using this code as it
        currently is.

        Returns:

        '''
        for param in self.parameters.keys():
            if 'Get' in self.parameters[param]:
                if param not in ['HSSpeed', 'VSSpeed', 'PreAmpGain', 'BitDepth', 'NumHSSpeed']:
                    self.GetParameter(param)
        if self.parameters['NumADChannels']['value'] == 1:
            self.GetParameter('BitDepth', 0)
            if 'value' in self.parameters['OutAmp']:
                self.GetParameter('NumHSSpeed', 0, self.parameters['OutAmp']['value'])

        return self.parameters

    # def SetAllParameters(self, parameters):
    #     # msg = ''
    #     for param in parameters.keys():
    #         if param not in ['EMMode', 'FastKinetics', 'ImageRotate', 'HSSpeed', 'VSSpeed', 'PreAmpGain', 'BitDepth',
    #                          'NumHSSpeed']:
    #             if 'Set' in self.parameters[param] and 'value' in self.parameters[param]:
    #                 if self.parameters[param]['value'] is not None:
    #                     if hasattr(self.parameters[param]['value'], '__iter__'):
    #                         self.SetParameter(param, *self.parameters[param]['value'])
    #                     else:
    #                         self.SetParameter(param, self.parameters[param]['value'])
    #
    #     #             else:
    #     #                 msg += param + ' '
    #     # self._logger.info('')

    '''Used functions'''

    def abort(self):
        try:
            self._dllWrapper('AbortAcquisition')
        except AndorWarning:
            pass

    def Initialize(self):
        self._dllWrapper('Initialize', outputs=(c_char(),))
        self.GetAllParameters()

        self.SetParameter('ReadMode', 4)
        self.SetParameter('AcquisitionMode', 1)
        self.SetParameter('TriggerMode', 0)
        self.SetParameter('Exposure', 0.01)
        self.SetParameter('Image', 1, 1, 1, self.parameters['DetectorShape']['value'][0], 1,
                          self.parameters['DetectorShape']['value'][1])
        self.SetParameter('Shutter', 1, 0, 1, 1)
        self.SetParameter('SetTemperature', -60)
        self.SetParameter('CoolerMode', 0)
        self.SetParameter('FanMode', 0)
        self.SetParameter('OutAmp', 0)

    # @background_action
    @locked_action
    def capture(self):
        """Capture function for Andor

        Wraps the three steps required for a camera acquisition: StartAcquisition, WaitForAcquisition and
        GetAcquiredData. The function also takes care of ensuring that the correct shape of array is passed to the
        GetAcquiredData call, according to the currently set parameters of the camera.

        Returns
        -------
        A numpy array containing the captured image(s)
        The number of images taken
        The shape of the images taken

        """
        self._dllWrapper('StartAcquisition')
        self._dllWrapper('WaitForAcquisition')
        self.WaitForDriver()

        if self.parameters['AcquisitionMode']['value'] == 4:
            num_of_images = 1  # self.parameters['FastKinetics']['value'][1]
            image_shape = (self.parameters['FastKinetics']['value'][-1], self.parameters['DetectorShape']['value'][0])
        else:
            if self.parameters['AcquisitionMode']['value'] == 1:
                num_of_images = 1
            elif self.parameters['AcquisitionMode']['value'] == 2:
                num_of_images = 1
            elif self.parameters['AcquisitionMode']['value'] == 3:
                num_of_images = self.parameters['NKin']['value']
            else:
                raise NotImplementedError('Acquisition Mode %g' % self.parameters['AcquisitionMode']['value'])

            if self.parameters['ReadMode']['value'] == 0:
                if self.parameters['IsolatedCropMode']['value'][0]:
                    image_shape = (
                    self.parameters['IsolatedCropMode']['value'][2] / self.parameters['IsolatedCropMode']['value'][4],)
                else:
                    image_shape = (self.parameters['DetectorShape']['value'][0] / self.parameters['FVBHBin']['value'],)
            elif self.parameters['ReadMode']['value'] == 3:
                image_shape = (self.parameters['DetectorShape']['value'][0],)
            elif self.parameters['ReadMode']['value'] == 4:
                if self.parameters['IsolatedCropMode']['value'][0]:
                    image_shape = (
                    self.parameters['IsolatedCropMode']['value'][1] / self.parameters['IsolatedCropMode']['value'][3],
                    self.parameters['IsolatedCropMode']['value'][2] / self.parameters['IsolatedCropMode']['value'][4])
                else:
                    image_shape = (
                        (self.parameters['Image']['value'][5] - self.parameters['Image']['value'][4] + 1) /
                        self.parameters['Image']['value'][1],
                        (self.parameters['Image']['value'][3] - self.parameters['Image']['value'][2] + 1) /
                        self.parameters['Image']['value'][0],)
            else:
                raise NotImplementedError('Read Mode %g' % self.parameters['ReadMode']['value'])

        dim = num_of_images * np.prod(image_shape)
        cimageArray = c_int * dim
        cimage = cimageArray()
        if '_logger' in self.__dict__:
            self._logger.debug('Getting AcquiredData for %i images with dimension %s' % (num_of_images, image_shape))
        try:
            self._dllWrapper('GetAcquiredData', inputs=({'type': c_int, 'value': dim},), outputs=(cimage,),
                             reverse=True)
            imageArray = []
            for i in range(len(cimage)):
                imageArray.append(cimage[i])
        except RuntimeWarning as e:
            if '_logger' in self.__dict__:
                self._logger.warn('Had a RuntimeWarning: %s' % e)
            imageArray = []
            for i in range(len(cimage)):
                imageArray.append(0)

        return imageArray, num_of_images, image_shape

    # @locked_action
    def SetImage(self, *params):
        """Set camera parameters for either the IsolatedCrop mode or Image mode

        Parameters
        ----------
        params  optional, inputs for either the IsolatedCrop mode or Image mode

        Returns
        -------

        """
        if self.parameters['IsolatedCropMode']['value'][0]:
            if len(params) == 0:
                params += (self.parameters['IsolatedCropMode']['value'])
            elif len(params) != 5:
                raise ValueError('Wrong number of parameters (need bool, cropheight, cropwidth, vbin, hbin')

            # Making sure we pass a valid set of parameters
            params = list(params)
            params[1] -= (params[1]) % params[3]
            params[2] -= (params[2]) % params[4]
            self.SetParameter('IsolatedCropMode', *params)
        else:
            if len(params) == 0:
                params = self.parameters['Image']['value']
            elif len(params) != 6:
                raise ValueError('Wrong number of parameters (need hbin, vbin, hstart, hend, vstart, vend')

            # Making sure we pass a valid set of parameters
            params = list(params)
            params[3] -= (params[3] - params[2] + 1) % params[0]
            params[5] -= (params[5] - params[4] + 1) % params[1]
            self.SetParameter('Image', *params)

    # @locked_action
    def SetROI(self, *params):
        if len(params) != 4:
            raise ValueError('Wrong number of inputs')
        current_binning = self.parameters['Image']['value'][:2]
        self.SetImage(*(current_binning + params))

    # @locked_action
    def SetBinning(self, *params):
        if len(params) not in [1, 2]:
            raise ValueError('Wrong number of inputs')
        if len(params) == 1:
            params += params
        current_ROI = self.parameters['Image']['value'][2:]
        self.SetImage(*(params + current_ROI))

    @locked_action
    def SetFastKinetics(self, n_rows=None):
        """Set the parameters for the Fast Kinetic mode

        Uses the already set parameters of exposure time, ReadMode, and binning as defaults to be passed to the Fast
        Kinetic parameter setter

        Parameters
        ----------
        n_rows

        Returns
        -------

        """

        if n_rows is None:
            n_rows = self.parameters['FastKinetics']['value'][0]

        series_Length = int(self.parameters['DetectorShape']['value'][1] / n_rows) - 1
        expT = self.parameters['AcquisitionTimings']['value'][0]
        mode = self.parameters['ReadMode']['value']
        hbin = self.parameters['Image']['value'][0]
        vbin = self.parameters['Image']['value'][1]
        offset = self.parameters['DetectorShape']['value'][1] - n_rows

        self.SetParameter('FastKinetics', n_rows, series_Length, expT, mode, hbin, vbin, offset)

    def GetStatus(self):
        error = self._dllWrapper('GetStatus', outputs=(c_int(),))
        self._status = ERROR_CODE[error]
        return self._status

    @locked_action
    def WaitForDriver(self):
        """
        This function is here because the dll.WaitForAcquisition does not work when in Accumulate mode

        Returns
        -------

        """
        status = c_int()
        self.dll.GetStatus(byref(status))
        while ERROR_CODE[status.value] == 'DRV_ACQUIRING':
            time.sleep(0.1)
            self.dll.GetStatus(byref(status))

    def CoolerON(self):
        self._dllWrapper('CoolerON')

    def CoolerOFF(self):
        self._dllWrapper('CoolerOFF')

    def IsCoolerOn(self):
        self.Cooler = self._dllWrapper('IsCoolerOn', outputs=(c_int(),))
        return self.Cooler

    def GetSeriesProgress(self):
        acc = c_long()
        series = c_long()
        error = self.dll.GetAcquisitionProgress(byref(acc), byref(series))
        if ERROR_CODE[error] == "DRV_SUCCESS":
            return series.value
        else:
            return None

    def GetAccumulationProgress(self):
        acc = c_long()
        series = c_long()
        error = self.dll.GetAcquisitionProgress(byref(acc), byref(series))
        if ERROR_CODE[error] == "DRV_SUCCESS":
            return acc.value
        else:
            return None


class Andor(Camera, AndorBase):
    metadata_property_names = ('Exposure', 'AcquisitionMode', 'TriggerMode')

    # def get_metadata(self, name):
    #     return self.GetParameter(name)
    #
    # def set_metadata(self, name, value):
    #     return self.SetParameter(name, value)

    def __init__(self, **kwargs):
        # Camera.__init__(self)
        AndorBase.__init__(self)
        Camera.__init__(self)
        # super(Andor, self).__init__()

        # self.wvl_to_pxl = kwargs['wvl_to_pxl']
        # self.magnification = kwargs['magnification']
        # self.pxl_size = kwargs['pxl_size']

        self.CurImage = None
        self.BGImage = None
        self.SubtractBG = None

        # self.acquiring = False
        self.isAborted = False
        # self.pause = False

    '''Used functions'''

    def Abort(self):
        self.isAborted = True
        self.abort()

    def raw_snapshot(self):
        try:
            imageArray, num_of_images, image_shape = self.capture()

            self.imageArray = imageArray

            # The image is reversed depending on whether you read in the conventional CCD register or the EM register, so we reverse it back
            if self.parameters['OutAmp']['value']:
                self.CurImage = np.reshape(self.imageArray, (num_of_images,) + image_shape)[..., ::-1]
            else:
                self.CurImage = np.reshape(self.imageArray, (num_of_images,) + image_shape)

            return self.CurImage[0], 1
        except Exception as e:
            self._logger.warn("Couldn't Capture because %s" % e)

    def get_camera_parameter(self, parameter_name):
        self.GetParameter(parameter_name)
    def set_camera_parameter(self, parameter_name, parameter_value):
        self.SetParameter(parameter_name, parameter_value)

    def get_qt_ui(self):
        return AndorUI(self)
    #
    # def getRelevantParameters(self):
    #     relevant_parameters = ['AcquisitionMode', 'ReadMode', 'Image', 'Exposure', 'NKin']
    #     dicc = {}
    #     for param in relevant_parameters:
    #         dicc[param] = self.parameters[param]['value']
    #     return dicc
    # def setRelevantParameters(self, parameters):
    #     for param in parameters:
    #         if hasattr(parameters[param], '__iter__'):
    #             self.SetParameter(param, *parameters[param])
    #         else:
    #             self.SetParameter(param, parameters[param])

    '''Not-used functions'''

    # def SaveAsBmp(self, path):
    #     im = Image.new("RGB", (self.parameters['DetectorShape'][0], self.parameters['DetectorShape'][1]), "white")
    #     pix = im.load()
    #
    #     for i in range(len(self.imageArray)):
    #         (row, col) = divmod(i, self.parameters['DetectorShape'][0])
    #         picvalue = int(round(self.imageArray[i] * 255.0 / 65535))
    #         pix[col, row] = (picvalue, picvalue, picvalue)
    #
    #     im.save(path, "BMP")
    #
    # def SaveAsTxt(self, path):
    #     file = open(path, 'w')
    #
    #     for line in self.imageArray:
    #         file.write("%g\n" % line)
    #
    #     file.close()
    #
    # def SaveAsBmpNormalised(self, path):
    #     im = Image.new("RGB", (self.parameters['DetectorShape'][0], self.parameters['DetectorShape'][1]), "white")
    #     pix = im.load()
    #
    #     maxIntensity = max(self.imageArray)
    #
    #     for i in range(len(self.imageArray)):
    #         (row, col) = divmod(i, self.parameters['DetectorShape'][0])
    #         picvalue = int(round(self.imageArray[i] * 255.0 / maxIntensity))
    #         pix[col, row] = (picvalue, picvalue, picvalue)
    #
    #     im.save(path, "BMP")
    #
    # def SaveAsFITS(self, filename, type):
    #     error = self.dll.SaveAsFITS(filename, type)
    #     self.verbose(ERROR_CODE[error], sys._getframe().f_code.co_name)
    #     return ERROR_CODE[error]


class AndorUI(QtWidgets.QWidget):
    def __init__(self, andor):
        assert isinstance(andor, Andor), "instrument must be an Andor"
        super(AndorUI, self).__init__()
        self.ImageUpdated = QtCore.SIGNAL('AndorImageUpdated')
        self.captureThread = None
        self.Andor = andor
        self.DisplayWidget = None

        uic.loadUi((os.path.dirname(__file__) + '/andor.ui'), self)

        self._setup_signals()
        self.updateGUI()
        self.BinningChanged()

        # self.Andor.updateGUI.connect(self.updateGUI)

    def __del__(self):
        self._stopTemperatureThread = True

    def _setup_signals(self):
        self.comboBoxAcqMode.activated.connect(self.AcquisitionModeChanged)
        self.comboBoxBinning.activated.connect(self.BinningChanged)
        self.comboBoxReadMode.activated.connect(self.ReadModeChanged)
        self.comboBoxTrigMode.activated.connect(self.TrigChanged)
        self.spinBoxNumFrames.valueChanged.connect(self.NumFramesChanged)
        self.spinBoxNumAccum.valueChanged.connect(self.NumAccumChanged)
        self.spinBoxNumRows.valueChanged.connect(self.NumRowsChanged)
        self.checkBoxROI.stateChanged.connect(self.ROI)
        self.checkBoxCrop.stateChanged.connect(self.IsolatedCrop)
        self.checkBoxCooler.stateChanged.connect(self.Cooler)
        # self.checkBoxAutoExp.stateChanged.connect(self.AutoExpose)
        self.checkBoxEMMode.stateChanged.connect(self.OutputAmplifierChanged)
        self.spinBoxEMGain.valueChanged.connect(self.EMGainChanged)
        self.lineEditExpT.returnPressed.connect(self.ExposureChanged)
        self.pushButtonDiv5.clicked.connect(lambda: self.ExposureChanged('/'))
        self.pushButtonTimes5.clicked.connect(lambda: self.ExposureChanged('x'))

        self.pushButtonCapture.clicked.connect(self.Capture)
        self.pushButtonLive.clicked.connect(self.Live)
        self.pushButtonAbort.clicked.connect(self.Abort)

    @background_action
    def _constantlyUpdateTemperature(self):
        self._stopTemperatureThread = False
        self.Andor.GetParameter('CurrentTemperature')
        while np.abs(self.Andor.parameters['CurrentTemperature']['value'] -
                             self.Andor.parameters['SetTemperature']['value']) > 2:
            temp = self.Andor.GetParameter('CurrentTemperature')
            self.checkBoxCooler.setText('Cooler (%g)' % temp)
            for ii in range(100):
                if self._stopTemperatureThread:
                    return
                time.sleep(0.1)

    # GUI FUNCTIONS
    def updateGUI(self):
        trig_modes = {0: 0, 1: 1, 6: 2}
        self.comboBoxAcqMode.setCurrentIndex(self.Andor.parameters['AcquisitionMode']['value'] - 1)
        self.comboBoxReadMode.setCurrentIndex(self.Andor.parameters['ReadMode']['value'])
        self.comboBoxTrigMode.setCurrentIndex(trig_modes[self.Andor.parameters['TriggerMode']['value']])
        self.comboBoxBinning.setCurrentIndex(np.log2(self.Andor.parameters['Image']['value'][0]))
        self.spinBoxNumFrames.setValue(self.Andor.parameters['NKin']['value'])

        self.Andor.GetParameter('AcquisitionTimings')
        self.lineEditExpT.setText(
            str(float('%#e' % self.Andor.parameters['AcquisitionTimings']['value'][0])).rstrip('0'))

    def Cooler(self):
        if self.checkBoxCooler.isChecked():
            if not self.Andor.IsCoolerOn():
                self.Andor.CoolerON()
                self.TemperatureUpdateThread = self._constantlyUpdateTemperature()
        else:
            if self.Andor.IsCoolerOn():
                self.Andor.CoolerOFF()
                if self.TemperatureUpdateThread.isAlive():
                    self._stopTemperatureThread = True

    def AcquisitionModeChanged(self):
        available_modes = ['Single', 'Accumulate', 'Kinetic', 'Fast Kinetic']
        currentMode = self.comboBoxAcqMode.currentText()
        self.Andor.SetParameter('AcquisitionMode', available_modes.index(currentMode) + 1)

    def ReadModeChanged(self):
        available_modes = ['FVB', 'Multi-track', 'Random track', 'Single track', 'Image']
        currentMode = self.comboBoxReadMode.currentText()
        self.Andor.SetParameter('ReadMode', available_modes.index(currentMode))

    def TrigChanged(self):
        available_modes = {'Internal': 0, 'External': 1, 'ExternalStart': 6}
        currentMode = self.comboBoxTrigMode.currentText()
        self.Andor.SetParameter('TriggerMode', available_modes[currentMode])

    def OutputAmplifierChanged(self):
        if self.checkBoxEMMode.isChecked():
            self.Andor.SetParameter('OutAmp', 0)
        else:
            self.Andor.SetParameter('OutAmp', 1)
        self.ROI()

    def BinningChanged(self):
        current_binning = int(self.comboBoxBinning.currentText()[0])
        if self.Andor.parameters['IsolatedCropMode']['value'][0]:
            params = list(self.Andor.parameters['IsolatedCropMode']['value'])
            params[3] = current_binning
            params[4] = current_binning
            print 'BinningChanged: ', params
            self.Andor.SetImage(*params)
        else:
            self.Andor.SetImage(current_binning, current_binning, *self.Andor.parameters['Image']['value'][2:])
        self.Andor.SetParameter('FVBHBin', current_binning)

    def NumFramesChanged(self):
        num_frames = self.spinBoxNumFrames.value()
        self.Andor.SetParameter('NKin', num_frames)

    def NumAccumChanged(self):
        num_frames = self.spinBoxNumAccum.value()
        # self.Andor.SetNumberAccumulations(num_frames)
        self.Andor.SetParameter('NAccum', num_frames)

    def NumRowsChanged(self):
        num_rows = self.spinBoxNumRows.value()
        if self.Andor.parameters['AcquisitionMode']['value'] == 4:
            self.Andor.SetFastKinetics(num_rows)
        elif self.Andor.parameters['ReadMode']['value'] == 3:
            self.Andor.SetParameter('SingleTrack', self.Andor.parameters['DetectorShape']['value'][1] / 2, num_rows)
        else:
            self.Andor._logger.info('Changing the rows only works in Fast Kinetic or in Single Track mode')

    def ExposureChanged(self, input=None):
        if input is None:
            expT = float(self.lineEditExpT.text())
        elif input == 'x':
            expT = float(self.lineEditExpT.text()) * 5
        elif input == '/':
            expT = float(self.lineEditExpT.text()) / 5

        # self.Andor.SetExposureTime(expT)
        self.Andor.SetParameter('Exposure', expT)
        # self.Andor.GetAcquisitionTimings()
        display_str = str(float('%#e' % self.Andor.parameters['AcquisitionTimings']['value'][0])).rstrip('0')
        self.lineEditExpT.setText(display_str)

    def EMGainChanged(self):
        gain = self.spinBoxEMGain.value()
        self.Andor.SetParameter('EMGain', gain)

    # @locked_action
    def IsolatedCrop(self):
        if self.DisplayWidget is None:
            return
        if hasattr(self.DisplayWidget, 'CrossHair1') and hasattr(self.DisplayWidget, 'CrossHair2'):
            current_binning = int(self.comboBoxBinning.currentText()[0])
            pos1 = self.DisplayWidget.CrossHair1.pos()
            pos2 = self.DisplayWidget.CrossHair2.pos()
            shape = self.Andor.parameters['DetectorShape']['value']
            maxx, minx = map(lambda x: int(x),
                             (min(pos1[0], pos2[0]), max(pos1[0], pos2[0])))
            miny, maxy = map(lambda x: int(x),  # shape[1] -
                             (min(pos1[1], pos2[1]), max(pos1[1], pos2[1])))
            if self.checkBoxCrop.isChecked():
                if self.checkBoxROI.isChecked():
                    self.checkBoxROI.setChecked(False)
                self.Andor.parameters['IsolatedCropMode']['value'] = (1,)
                self.Andor.SetImage(1, maxy, minx, current_binning, current_binning)
                # self.Andor.SetParameter('IsolatedCropMode', 1, maxy, minx, current_binning, current_binning)
            else:
                self.Andor.SetParameter('IsolatedCropMode', 0, maxy, minx, current_binning, current_binning)
                self.Andor.SetImage()

    # @locked_action
    def ROI(self):
        if self.DisplayWidget is None:
            return
        if hasattr(self.DisplayWidget, 'CrossHair1') and hasattr(self.DisplayWidget, 'CrossHair2'):
            hbin, vbin = self.Andor.parameters['Image']['value'][:2]
            if self.checkBoxROI.isChecked():
                if self.checkBoxCrop.isChecked():
                    self.checkBoxCrop.setChecked(False)
                pos1 = self.DisplayWidget.CrossHair1.pos()
                pos2 = self.DisplayWidget.CrossHair2.pos()
                shape = self.Andor.parameters['DetectorShape']['value']
                # print 'GUI ROI. CrossHair: ', pos1, pos2
                maxx, minx = map(lambda x: shape[0] - int(x),
                                 (min(pos1[0], pos2[0]), max(pos1[0], pos2[0])))
                miny, maxy = map(lambda x: int(x) + 1,  # shape[1] -
                                 (min(pos1[1], pos2[1]), max(pos1[1], pos2[1])))

                # print 'ROI. ImageInfo: ',hbin, vbin, minx, maxx, miny, maxy
                # if self.Andor.parameters['OutAmp']['value']:
                #     self.Andor.SetImage(hbin, vbin, shape[0]-maxx, shape[0]-minx, miny, maxy)
                # else:
                self.Andor.SetImage(hbin, vbin, minx, maxx, miny, maxy)
                # self.Andor.parameters['Image']['value'] = [hbin, vbin, minx, miny, maxx, maxy]
                # print 'GUI ROI. ShapeInfo: ', self.Andor.parameters['Image']
                # self.Andor.SetImage()
            else:
                self.Andor.SetParameter('Image', hbin, vbin, 1, self.Andor.parameters['DetectorShape']['value'][0],
                                        1, self.Andor.parameters['DetectorShape']['value'][1])
                # self.Andor.parameters['Image'] = [1, 1, self.Andor.parameters['DetectorShape'][0],
                #                                   self.Andor.parameters['DetectorShape'][1]]
                # self.Andor.SetImage()

    def Capture(self, wait=True):
        # print 'Capture'
        # t = self.Andor.Capture()
        if self.captureThread is not None:
            if not self.captureThread.isFinished():
                return
        self.captureThread = CaptureThread(self.Andor)
        self.connect(self.captureThread, self.captureThread.updateImage, self.updateImage)
        # self.captureThread.finished.connect(self.updateImage)
        self.captureThread.start()

        if wait:
            self.captureThread.wait()

    def Live(self, wait=True):
        if self.captureThread is not None:
            if not self.captureThread.isFinished():
                return
        self.captureThread = CaptureThread(self.Andor, live=True)
        self.connect(self.captureThread, self.captureThread.updateImage, self.updateImage)
        # self.captureThread.finished.connect(self.updateImage)
        self.captureThread.start()

        if wait:
            self.captureThread.wait()

    def Abort(self):
        # self.Andor.abort = True
        self.Andor.Abort()

    def updateImage(self):
        if self.DisplayWidget is None:
            from Experiments import maingui
            self.DisplayWidget = maingui.DisplayWidget()
        if self.DisplayWidget.isHidden():
            self.DisplayWidget.show()

        if len(self.Andor.CurImage.shape) == 2:
            self.DisplayWidget.splitter.setSizes([0, 1])
            if self.Andor.CurImage.shape[0] > self.DisplayWidget._max_num_line_plots:
                self.Andor._logger.warn('Trying to display too many lines')
                number = self.DisplayWidget._max_num_line_plots
            else:
                number = self.Andor.CurImage.shape[0]
            for ii in range(number):
                self.DisplayWidget.plot[ii].setData(self.Andor.CurImage[ii])
        else:
            self.DisplayWidget.splitter.setSizes([1, 0])
            image = np.transpose(self.Andor.CurImage, (0, 2, 1))
            offset = ((self.Andor.parameters['DetectorShape']['value'][0] - self.Andor.parameters['Image']['value'][3]),
                      self.Andor.parameters['Image']['value'][4] - 1)
            if self.Andor.parameters['IsolatedCropMode']['value'][0]:
                scale = self.Andor.parameters['IsolatedCropMode']['value'][-2:]
            else:
                scale = self.Andor.parameters['Image']['value'][:2]
            if image.shape[0] == 1:
                image = image[0]
                self.DisplayWidget.ImageDisplay.setImage(image,
                                                         pos=offset, autoRange=False,
                                                         scale=scale)
            else:
                self.DisplayWidget.ImageDisplay.setImage(image, xvals=0.99 * np.linspace(0, image.shape[0] - 1,
                                                                                         image.shape[0]),
                                                         pos=offset, autoRange=False,
                                                         scale=scale)
        self.emit(self.ImageUpdated, 'Andor')


class CaptureThread(QtCore.QThread):
    def __init__(self, andor, live=False):
        QtCore.QThread.__init__(self, parent=None)
        self.updateImage = QtCore.SIGNAL("UpdateImage")
        self.Andor = andor
        self.live = live

    def stop(self):
        self.isAborted = True
        self.wait()

    def run(self):
        if self.live:
            self.Andor.isAborted = False
            while not self.Andor.isAborted:
                try:
                    self.SingleAcquire()
                except AndorWarning:
                    pass
        else:
            self.SingleAcquire()
        self.Andor.isAborted = False

    def SingleAcquire(self):
        self.Andor.raw_snapshot()
        if self.Andor.parameters['AcquisitionMode']['value'] in [1, 2] and self.Andor.parameters['NKin']['value'] > 1:
            if self.Andor.parameters['SoftwareWaitBetweenCaptures']:
                time.sleep(self.Andor.parameters['SoftwareWaitBetweenCaptures'])

            final_array = np.zeros(
                (self.Andor.parameters['NKin']['value'],) + self.Andor.CurImage.shape[1:])
            final_array[0] = self.Andor.CurImage[0]
            for ii in range(1, self.Andor.parameters['NKin']['value']):
                if self.Andor.isAborted:
                    break
                self.Andor.Capture()
                final_array[ii] = self.Andor.CurImage[0]
            self.Andor.CurImage = final_array

        self.emit(self.updateImage)


class WaitThread(QtCore.QThread):
    def __init__(self, andor):
        QtCore.QThread.__init__(self, parent=None)
        self.Andor = andor

    def run(self):
        self.Andor._logger.infor('Waiting for temperature to come up')
        temp = 30
        try:
            temp = self.Andor._dllWrapper('GetTemperature', outputs=(c_int(),))[0]
        except AndorWarning as warn:
            if warn.error_name != 'DRV_TEMP_OFF':
                raise warn
        if self.Andor.IsCoolerOn():
            self.Andor.CoolerOFF()
        if temp < 30:
            toggle = windll.user32.MessageBoxA(0, 'Camera is cold (%g), do you want to wait before ShutDown? '
                                                  '\n Not waiting can cause irreversible damage' % temp, '', 4)
            if toggle == 7:
                return
            else:
                while temp < -20:
                    self.Andor._logger.info('Waiting for temperature to come up. %g' % temp)
                    time.sleep(10)
                    try:
                        temp = self.Andor._dllWrapper('GetTemperature', outputs=(c_int(),))[0]
                    except AndorWarning as warn:
                        if warn.error_name != 'DRV_TEMP_OFF':
                            raise warn


ERROR_CODE = {
    20001: "DRV_ERROR_CODES",
    20002: "DRV_SUCCESS",
    20003: "DRV_VXNOTINSTALLED",
    20006: "DRV_ERROR_FILELOAD",
    20007: "DRV_ERROR_VXD_INIT",
    20010: "DRV_ERROR_PAGELOCK",
    20011: "DRV_ERROR_PAGE_UNLOCK",
    20013: "DRV_ERROR_ACK",
    20024: "DRV_NO_NEW_DATA",
    20026: "DRV_SPOOLERROR",
    20034: "DRV_TEMP_OFF",
    20035: "DRV_TEMP_NOT_STABILIZED",
    20036: "DRV_TEMP_STABILIZED",
    20037: "DRV_TEMP_NOT_REACHED",
    20038: "DRV_TEMP_OUT_RANGE",
    20039: "DRV_TEMP_NOT_SUPPORTED",
    20040: "DRV_TEMP_DRIFT",
    20050: "DRV_COF_NOTLOADED",
    20053: "DRV_FLEXERROR",
    20066: "DRV_P1INVALID",
    20067: "DRV_P2INVALID",
    20068: "DRV_P3INVALID",
    20069: "DRV_P4INVALID",
    20070: "DRV_INIERROR",
    20071: "DRV_COERROR",
    20072: "DRV_ACQUIRING",
    20073: "DRV_IDLE",
    20074: "DRV_TEMPCYCLE",
    20075: "DRV_NOT_INITIALIZED",
    20076: "DRV_P5INVALID",
    20077: "DRV_P6INVALID",
    20083: "P7_INVALID",
    20089: "DRV_USBERROR",
    20091: "DRV_NOT_SUPPORTED",
    20095: "DRV_INVALID_TRIGGER_MODE",
    20099: "DRV_BINNING_ERROR",
    20990: "DRV_NOCAMERA",
    20991: "DRV_NOT_SUPPORTED",
    20992: "DRV_NOT_AVAILABLE"
}


def test_all():
    import matplotlib.pyplot as plt
    andor = Andor() #wvl_to_pxl=32.5 / 1600, magnification=30, pxl_size=16)
    print andor
    andor.show_gui(True)
    # params = andor.GetAllParameters()
    # andor.SetAllParameters(params)
    print 'Done'

    # data = []
    # andor.StartAcquisition()
    # andor.GetAcquiredData(data)
    #
    # # print data
    # print type(data)
    # print type(andor.CurImage)
    # print andor.CurImage[0:10]
    #
    # plt.imshow(andor.CurImage)
    # plt.show()


def test_GUI():
    from gui import QtWidgets
    import sys

    andor = Andor(verbosity=True, idn='AndorDebug', wvl_to_pxl=32.5 / 1600, magnification=30, pxl_size=16)

    app = QtWidgets.QApplication(sys.argv)
    app.aboutToQuit.connect(app.deleteLater)

    # # send the experiment instance to the main dialog
    dialog = AndorUI(andor)
    dialog.show()

    try:
        from IPython.lib.guisupport import start_event_loop_qt4
        start_event_loop_qt4(app)
    except ImportError:
        app.exec_()

    sys.exit(app.exec_())


if __name__ == '__main__':
    test_all()