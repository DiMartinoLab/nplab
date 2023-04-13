"""
Issues:
    - The waitStop property for moving doesn't really work, and if you send two move commands quickly after each other,
    the system doesn't react fast enough and doesn't reach the final destination.
    
Modified by Trung 16.02.23 for rig2 SM100 controller and Newport stages
"""
from __future__ import print_function

from builtins import map
from builtins import hex
from builtins import str
from builtins import range
import time

from nplab.instrument.serial_instrument import SerialInstrument
from nplab.instrument.stage import Stage

# never wait for more than this e.g. during wait_states
MAX_WAIT_TIME_SEC = 300

# time to wait after sending a command. This number has been arrived at by
# trial and error
COMMAND_WAIT_TIME_SEC = 1

# States from page 65 of the manual
STATE_NOT_REFERENCED_FROM_RESET = '0A'
STATE_NOT_REFERENCED_FROM_CONFIGURATION = '0C'
STATE_READY_FROM_HOMING = '32'
STATE_READY_FROM_MOVING = '33'

STATE_CONFIGURATION = '14'

STATE_DISABLE_FROM_READY = '3C'
STATE_DISABLE_FROM_MOVING = '3D'
STATE_DISABLE_FROM_JOGGING = '3E'


class SMC100ReadTimeOutException(Exception):
    def __init__(self):
        super(SMC100ReadTimeOutException, self).__init__('Read timed out')


class SMC100WaitTimedOutException(Exception):
    def __init__(self):
        super(SMC100WaitTimedOutException, self).__init__('Wait timed out')


class SMC100DisabledStateException(Exception):
    def __init__(self, state):
        super(SMC100DisabledStateException, self).__init__('Disabled state encountered: ' + state)


class SMC100RS232CorruptionException(Exception):
    def __init__(self, c):
        super(SMC100RS232CorruptionException, self).__init__('RS232 corruption detected: %s' % (hex(ord(c))))


class SMC100InvalidResponseException(Exception):
    def __init__(self, cmd, resp):
        s = 'Invalid response to %s: %s' % (cmd, resp)
        super(SMC100InvalidResponseException, self).__init__(s)

        
        
class SMC100(SerialInstrument):
    """
    Class to interface with Newport's SMC100 controller.
    The SMC100 accepts commands in the form of:
      <ID><command><arguments><CR><LF>
    Reply, if any, will be in the form
      <ID><command><result><CR><LF>
    There is minimal support for manually setting stage parameter as Newport's
    ESP stages can supply the SMC100 with the correct configuration parameters.
    Some effort is made to take up backlash, but this should not be trusted too
    much.
    The move commands must be used with care, because they make assumptions
    about the units which is dependent on the STAGE. I only have TRB25CC, which
    has native units of mm. A more general implementation will move the move
    methods into a stage class.
    """

    def __init__(self, port, smcID=(1, ), **kwargs):
        """
        If backlash_compensation is False, no backlash compensation will be done.
        If silent is False, then additional output will be emitted to aid in
        debugging.
        If sleepfunc is not None, then it will be used instead of time.sleep. It
        will be given the number of seconds (float) to sleep for, and is provided
        for ease integration with single threaded GUIs.
        Note that this method only connects to the controller, it otherwise makes
        no attempt to home or configure the controller for the attached stage. This
        delibrate to minimise realworld side effects.
        If the controller has previously been configured, it will suffice to simply
        call home() to take the controller out of not referenced mode. For a brand
        new controller, call reset_and_configure().
        """
        self.port_settings = dict(baudrate=57600,
                    bytesize=8,
                    stopbits=1,
                    parity='N',
                    xonxoff=True,
                    timeout=0.050)

        SerialInstrument.__init__(self, port)
        

        # self._logger.debug('Connecting to SMC100 on %s' % (port))

        self.software_home = None
        self._last_sendcmd_time = 0
        if not hasattr(smcID, '__iter__'):
            smcID = (smcID, )
        self._smcID = list(smcID)
        self.axis_names = ()
        print('found controllers: ')
        for id in self._smcID:
            self.axis_names += (str(id), )
            self._send_cmd('ID', id, '?', True)  # Just testing the connection
            
            

    def __del__(self):
        self.close()

    def _send_cmd(self, command, axes=None, argument=None, expect_response=False, retry=False):
        """
        Send the specified command along with the argument, if any. The response
        is checked to ensure it has the correct prefix, and is returned WITHOUT
        the prefix.
        It is important that for GET commands, e.g. 1ID?, the ? is specified as an
        ARGUMENT, not as part of the command. Doing so will result in assertion
        failure.
        If expect_response is True, a response is expected from the controller
        which will be verified and returned without the prefix.
        If expect_response is True, and retry is True or an integer, then when the
        response does not pass verification, the command will be sent again for
        retry number of times, or until success if retry is True.
        The retry option MUST BE USED CAREFULLY. It should ONLY be used read-only
        commands, because otherwise REPEATED MOTION MIGHT RESULT. In fact some
        commands are EXPLICITLY REJECTED to prevent this, such as relative move.
        """
        if axes is None:
            axes = self.axis_names #self._smcID[0]
        elif not hasattr(axes, '__iter__'):
            axes = (axes, )

        reply = ()
        for axis in axes:
            if type(axis) != str:
                axis = str(axis)
            assert command[-1] != '?'

            if argument is None:
                argument = ''

            prefix = axis + command
            tosend = prefix + str(argument)

            # prevent certain commands from being retried automatically
            no_retry_commands = ['PR', 'OR']
            if command in no_retry_commands:
                retry = False

            done = False
            while not done:
                if expect_response:
                    self.ser.flushInput()

                self.ser.flushOutput()

                self.ser.write(bytes(tosend,'utf-8'))
                self.ser.write(bytes('\r\n','utf-8')) #converted to bytes by sunny

                self.ser.flush()

                if expect_response:
                    try:
                        response = self._readline()
                        if response.startswith(prefix):
                            reply += (response[len(prefix):], )
                            done = True
                        else:
                            raise SMC100InvalidResponseException(command, response)
                    except Exception as ex:
                        if not retry or retry <= 0:
                            raise ex
                        else:
                            if type(retry) == int:
                                retry -= 1
                            continue
                else:
                    # we only need to delay when we are not waiting for a response
                    now = time.time()
                    dt = now - self._last_sendcmd_time
                    dt = COMMAND_WAIT_TIME_SEC - dt
                    # print dt
                    if dt > 0:
                        time.sleep(dt)
                    self._last_sendcmd_time = now
                    done = True
                    # return None
        #print(str(reply))
        return reply

    def _readline(self):
        """
        Returns a line, that is reads until \r\n.
        OK, so you are probably wondering why I wrote this. Why not just use
        self.ser.readline()?
        I am glad you asked.
        With python < 2.6, pySerial uses serial.FileLike, that provides a readline
        that accepts the max number of chars to read, and the end of line
        character.
        With python >= 2.6, pySerial uses io.RawIOBase, whose readline only
        accepts the max number of chars to read. io.RawIOBase does support the
        idea of a end of line character, but it is an attribute on the instance,
        which makes sense... except pySerial doesn't pass the newline= keyword
        argument along to the underlying class, and so you can't actually change
        it.
        """
        done = False
        line = str()
        # print 'reading line',
        while not done:
            c = self.ser.read().decode('utf-8') #decode added by sunny
            #print(c)
            # ignore \r since it is part of the line terminator
            if len(c) == 0:
                raise SMC100ReadTimeOutException()
            elif c == '\r':
                continue
            elif c == '\n':
                done = True
            elif ord(c) > 32 and ord(c) < 127:
                line += c
            else:
                raise SMC100RS232CorruptionException(c)

        return line

    def _wait_states(self, targetstates, ignore_disabled_states=False):
        """
        Waits for the controller to enter one of the the specified target state.
        Controller state is determined via the TS command.
        If ignore_disabled_states is True, disable states are ignored. The normal
        behaviour when encountering a disabled state when not looking for one is
        for an exception to be raised.
        Note that this method will ignore read timeouts and keep trying until the
        controller responds.  Because of this it can be used to determine when the
        controller is ready again after a command like PW0 which can take up to 10
        seconds to execute.
        If any disable state is encountered, the method will raise an error,
        UNLESS you were waiting for that state. This is because if we wait for
        READY_FROM_MOVING, and the stage gets stuck we transition into
        DISABLE_FROM_MOVING and then STAY THERE FOREVER.
        The state encountered is returned.
        """
        starttime = time.time()
        done = [False]*len(self.axis_names)
        self._logger.debug('waiting for states %s' % (str(targetstates)))
        while not all(done):
            for axes in range(len(self.axis_names)):
                waittime = time.time() - starttime
                if waittime > MAX_WAIT_TIME_SEC:
                    raise SMC100WaitTimedOutException()

                try:
                    state = self.get_status()[axes][1]
                    if state in targetstates:
                        self._logger.debug('in state %s' % (state))
                        done[axes] = True
                        # return state
                    elif not ignore_disabled_states:
                        disabledstates = [
                            STATE_DISABLE_FROM_READY,
                            STATE_DISABLE_FROM_JOGGING,
                            STATE_DISABLE_FROM_MOVING]
                        if state in disabledstates:
                            raise SMC100DisabledStateException(state)

                except SMC100ReadTimeOutException:
                    self._logger.info('Read timed out, retrying in 1 second')
                    time.sleep(1)
                    continue

    def reset_and_configure(self):
        """
        Configures the controller by resetting it and then asking it to load
        stage parameters from an ESP compatible stage. This is then followed
        by a homing action.
        """
        self._send_cmd('RS')
        self._send_cmd('RS')

        self._wait_states(STATE_NOT_REFERENCED_FROM_RESET, ignore_disabled_states=True)

        stage = str(self._send_cmd('ID', None,'?', True))
        self._logger.info('Found stage %s' %stage)
        print(str(stage))
        # enter config mode
        self._send_cmd('PW1',1) #'PW1' added by sunny
        self._send_cmd('PW1',2)
        self._send_cmd('PW1',3) #'PW3' #added by Asia for 3rd axis
        self._wait_states(STATE_CONFIGURATION)
        
        # load stage parameters
        self._send_cmd('ZX2', 1)
        time.sleep(2)
        # enable stage ID check
        self._send_cmd('ZX2', 2)
        time.sleep(2)
        
        self._send_cmd('ZX2', 3) #added by Asia for 3rd axis, seems not necessary
        #time.sleep(2) #added by Asia for 3rd axis, seems not necessary
        
        # for Rig2
        self._send_cmd('BA', 1, '0.00862') #modified by Trung 16.02.23 - backlash values, printed on actuators
        self._send_cmd('BA', 2, '0.00853') #modified by Trung 16.02.23  - backlash values, printed on actuators
        self._send_cmd('BA', 3, '0.01151') #added by Asia - backlash values, printed on actuators        

        
        # exit configuration mode
        self._send_cmd('PW0', 1)
        # xtn20: wait for communication with stage
        time.sleep(5)
        self._send_cmd('PW0', 2)
        time.sleep(5) #added by Asia for 3rd axis
        self._send_cmd('PW0', 3) #added by Asia for 3rd axis
        self._wait_states(STATE_NOT_REFERENCED_FROM_CONFIGURATION)
        #homing part to make box flash green.#Asia: homing and putting in ready state
        self._send_cmd('OR',1)
        self._send_cmd('OR',2)
        self._send_cmd('OR',3) #added by Asia for 3rd axis

    def get_position(self, axis=None):
        pos = self._send_cmd('TP', axes=axis, argument='?', expect_response=True, retry=10)
        pos = list(map(float, pos))
        return pos

    def home(self, **kwargs):
        """
        Homes the controller. If waitStop is True, then this method returns when
        homing is complete.
        Note that because calling home when the stage is already homed has no
        effect, and homing is generally expected to place the stage at the
        origin, an absolute move to 0 um is executed after homing. This ensures
        that the stage is at origin after calling this method.
        Calling this method is necessary to take the controller out of not referenced
        state after a restart.
        """
        self._send_cmd('OR')
        if 'waitStop' in kwargs and kwargs['waitStop']:
            # wait for the controller to be ready
            st = self._wait_states((STATE_READY_FROM_HOMING, STATE_READY_FROM_MOVING))
            if st == STATE_READY_FROM_MOVING:
                self.move([0]*len(self.axis_names), **kwargs)
        else:
            self.move([0]*len(self.axis_names), **kwargs)

    def stop(self):
        self._send_cmd('ST')

    def get_status(self):
        """
        Executes TS? and returns the the error code as integer and state as string
        as specified on pages 64 - 65 of the manual.
        """

        resps = self._send_cmd('TS', argument='?', expect_response=True, retry=10)
        reply = ()
        for resp in resps:
            errors = int(resp[0:4], 16)
            state = resp[4:]
            assert len(state) == 2
            reply += ([errors, state], )

        return reply

    def move(self, pos, axis=None, relative=False, waitStop=True):
        if axis is None:
            axis = self.axis_names
        if not hasattr(pos, '__iter__'):
            pos = [pos]
        if relative:
            index = 0
            for ax in axis:
                self._send_cmd('PR', axes=ax, argument=pos[index]*1000)
                index += 1
        else:
            index = 0
            for ax in axis:
                self._send_cmd('PA', axes=ax, argument=pos[index])
                index += 1

        if waitStop:
            # If we were previously homed, then something like PR0 will have no
            # effect and we end up waiting forever for ready from moving because
            # we never left ready from homing. This is why STATE_READY_FROM_HOMING
            # is included.
            self._wait_states((STATE_READY_FROM_MOVING, STATE_READY_FROM_HOMING))

    def move_referenced(self, position_mm, **kwargs):
        """
        Moves to an absolute position referenced from the software home

        Args:
            position_mm: position from the software home
            **kwargs: kwargs to be passed to the move command

        Returns:

        """

        if not hasattr(position_mm, '__iter__'):
            position_mm = (position_mm, )

        final_pos = list(map(lambda x, y: x+y, self.software_home, position_mm))

        self.move(final_pos, **kwargs)

    def set_software_home(self):
        """
        Sets a software home, so that we can easily go back to similar sample positions

        Returns:

        """
        self.software_home = self.get_position()

    def go_software_home(self):
        self.move_referenced([0]*len(self.axis_names))

    def set_velocity(self, velocity):
        self._send_cmd('VA_Set', velocity)


if __name__ == '__main__':
    smc100 = SMC100('COM1', (1,2,3))
    smc100.show_gui()











