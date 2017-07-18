__author__ = 'alansanders, chrisgrosse'

from nplab.instrument.visa_instrument import VisaInstrument
from nplab.instrument.stage import Stage, StageUI
import time
import numpy as np
from functools import partial



class PIStageError(Exception):
    def __init__(self, error_code):
        self.error_code = error_code

    def __str__(self):
        return "Code %s: %s" % (self.error_code, self.get_error_message(self.error_code))

    def get_error_message(self, error_code):
        print "ERROR CODE:", error_code
        if (error_code == '1'):
            msg = 'PARAMETER SYNTAX ERROR!'
        elif (error_code == '5'):
            msg = 'CANNOT SET POSITION WHEN SERVO IS OFF!'
        elif (error_code == '10'):
            msg = 'CONTROLLER WAS STOPPED!'
        elif (error_code == '18'):
            msg = 'INVALID MACRO NAME!'
        elif (error_code == '19'):
            msg = 'ERROR WHILE RECORDING MACRO!'
        elif (error_code == '20'):
            msg = 'MACRO NOT FOUND!'
        elif (error_code == '23'):
            msg = 'ILLEGAL AXIS IDENTIFIER!'
        elif (error_code == '26'):
            msg = 'PARAMETER MISSING!'
        elif (error_code == '301'):
            msg = 'SEND BUFFER OVERFLOW!'
        elif (error_code == '302'):
            msg = 'VOLTAGE OUT OF LIMITS!'
        elif (error_code == '303'):
            msg = 'CANNOT SET VOLTAGE WHEN SERVO-MODE ON!'
        elif (error_code == '304'):
            msg = 'RECEIVED COMMAND IN TOO LONG!'
        elif (error_code == '307'):
            msg = 'TIMEOUT WHILE RECEIVING COMMAND'
        elif (error_code == '309'):
            msg = 'INSUFFICIENT SPACE TO STORE MACRO'
        else:
            msg = 'UNKNOWN ERROR'
        return msg


class PIStage(VisaInstrument, Stage):
    """
    Control interface for PI stages.
    """
    def __init__(self, address='ASRL3::INSTR',timeout = 10,baud_rate = 57600):
        super(PIStage, self).__init__(address=address)
        self.instr.read_termination = '\n'
        self.instr.write_termination = '\n'
        self.instr.baud_rate = 57600
   #     self.instr.timeout = 10
        self.axis_names = ('a', 'b')
        self.positions = [0 for ch in xrange(3)]
        self._stage_id = None
        #self.startup()

    def move(self, pos, axis=None, relative=False):
        if relative:
            self.set_axis_param(partial(self.move_axis, relative=True), pos, axis)
        else:
            self.set_axis_param(self.move_axis, pos, axis)

    def move_axis(self, pos, axis, relative=False):
        if not self.get_loop_mode(axis): # check whether servo-control mode is on
            self.set_loop_mode(1,axis)
        if relative:
            self.write('mvr {0}{1}'.format(axis, 1e6*pos))
        else:
            self.write('mov {0}{1}'.format(axis, 1e6*pos))
        self.wait_until_stopped(axis)
        self.check_status()

    def stop(self):
        """
        stops any ongoing movement of the positioner
        """
        self.write('#24')

    def get_position(self, axis=None):
        txt = self.get_axis_param(lambda axis: 1e-6*float(self.float_query('pos? {0}'.format(axis))), axis)
        self.check_status()
        return txt
    position = property(fget=get_position, doc="Current position of the stage")

    def get_voltage(self, axis=None):
        txt = self.get_axis_param(lambda axis: self.float_query('vol? {0}'.format(axis)),axis)
        self.check_status()
        return txt
    def set_voltage(self, voltage, axis=None):
        if self.get_loop_mode(axis): # servo-control mode needs to be switched off to set voltage!
            self.set_loop_mode(0,axis)
        self.set_axis_param(lambda value, axis: self.write('sva {0}{1}'.format(axis, voltage)), voltage, axis)
        self.check_status()
    piezo_voltages = property(fget=get_voltage, fset=set_voltage, doc="Current voltage of the stage")

    def check_status(self):
        error_code = self.query('err?')
        if (error_code != '0'):
            raise PIStageError(error_code)
        return error_code


    def is_moving(self, axes=None):
        """
        Returns True if any of the specified axes are in motion.
        In this case the position is polled 3 times and see if the stage stays close to
        its initial position.
        """
        positions = np.zeros((3, len(axes)))
        for i in xrange(3):
            positions[i] = [self.get_position(axis) for axis in axes]
            time.sleep(0.005)
        sum_of_diffs = np.sum(positions-positions[0], axis=1)
        if np.any(sum_of_diffs > 0.01):
            print sum_of_diffs
            return True
        else:
            return False

    def wait_until_stopped(self, axes=None): # same as in super class Stage
        """Block until the stage is no longer moving."""
        while self.is_moving(axes=axes):
            time.sleep(0.01)

    def startup(self):
        self.online = 1
        while not self.online:
            print self.online
        self.loop_mode = 1
        self.speed_mode = 0
        self.velocity = 100
        self.drift_compensation = 0
        self.instr.write('cto 132')
        self.instr.write('cto 232')
        self.instr.write('cto 332')

    def shutdown(self):
        self.loop_mode = 0
        self.online = 0

    def get_min_voltage(self, axis=None):
        txt = self.get_axis_param(lambda axis: float(self.float_query('vmi? {0}'.format(axis))), axis)
        self.check_status()
        return txt
    def set_min_voltage(self, voltage, axis=None):
        self.set_axis_param(lambda voltage, axis: self.write('vmi {0}{1}'.format(axis, voltage)), voltage, axis)
        self.check_status()

    def get_max_voltage(self, axis=None):
        txt = self.get_axis_param(lambda axis: float(self.float_query('vma? {0}'.format(axis))), axis)
        self.check_status()
        return txt
    def set_max_voltage(self, voltage, axis=None):
        self.set_axis_param(lambda voltage, axis: self.write('vma {0}{1}'.format(axis, voltage)), voltage, axis)
        self.check_status()

    def get_velocity(self, axis=None):
        return self.get_axis_param(lambda axis: float(self.float_query('vel? {0}'.format(axis))), axis)

    def set_velocity(self, value, axis=None):
        self.set_axis_param(lambda value, axis: self.write('vel {0}{1}'.format(axis, value)), value, axis)
        self.check_status()
    velocity = property(get_velocity, set_velocity)

    def get_drift_compensation(self, axis=None):
        return self.get_axis_param(lambda axis: bool(self.int_query('dco? {0}'.format(axis))), axis)
    def set_drift_compensation(self, value, axis=None):
        self.set_axis_param(lambda value, axis: self.write('dco {0}{1}'.format(axis, value)), value, axis)
        self.check_status()
    drift_compensation = property(get_drift_compensation, set_drift_compensation)

    def get_loop_mode(self, axis=None):
        txt = self.get_axis_param(lambda axis: bool(self.int_query('svo? {0}'.format(axis))), axis)
        self.check_status()
        return txt
    def set_loop_mode(self, value, axis=None):
        """
        Set the mode of each axis control loop
        :param value: servo control mode - 1 for closed loop, 0 for open loop
        :param axis:
        :return:
        """
        self.set_axis_param(lambda value, axis: self.write('svo {0}{1}'.format(axis, value)), value, axis)
        self.check_status()
    loop_mode = property(get_loop_mode, set_loop_mode)

    def get_speed_mode(self, axis=None):
        return self.get_axis_param(lambda axis: bool(self.int_query('vco? {0}'.format(axis))), axis)
    def set_speed_mode(self, value, axis=None):
        """
        Set the mode of each axis control loop
        :param value: speed control mode - 1 for controlled speed, 0 for fastest
        :param axis:
        :return:
        """
        self.set_axis_param(lambda value, axis: self.write('vco {0}{1}'.format(axis, value)), value, axis)
        self.check_status()
    speed_mode = property(get_speed_mode, set_speed_mode)

    def get_online(self):
        return bool(self.int_query('onl?'))
    def set_online(self, value):
        self.write('onl {0}'.format(value))
        self.check_status()
    online = property(get_online, set_online)

    def get_on_target(self):
       return bool(self.int_query('ont?'))
    on_target = property(get_on_target)

    def get_id(self):
        if self._stage_id is None:
            self._stage_id = self.query('*idn?')
        return self._stage_id
    stage_id = property(get_id)

    def get_qt_ui(self):
        return StageUI(self, stage_step_min=0.1e-9, stage_step_max=100e-6)


if __name__ == '__main__':
    stage = PIStage()
  #  stage.move((5e-6, 10e-6))
#    print stage.position
#    print stage.get_position()
#    print stage.get_position(axis=('a', 'b'))
#
#    import sys
#    from nplab.utils.gui import get_qt_app
#    app = get_qt_app()
#    ui = stage.get_qt_ui()
#    ui.show()
#    sys.exit(app.exec_())