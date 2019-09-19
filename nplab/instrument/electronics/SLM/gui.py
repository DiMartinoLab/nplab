# -*- coding: utf-8 -*-

from nplab.utils.gui import QtWidgets, QtGui, QtCore, uic, get_qt_app
import numpy as np
import os


class gratingsUi(QtWidgets.QWidget):
    def __init__(self, slm_gui):
        super(gratingsUi, self).__init__()
        uic.loadUi(os.path.join(os.path.dirname(__file__), 'ui_gratings.ui'), self)
        self.slm_gui = slm_gui
        self._connect()

    def _connect(self):
        self.pushButton_center.clicked.connect(lambda: self.update_gratings('center'))
        self.pushButton_up.clicked.connect(lambda: self.update_gratings('up'))
        self.pushButton_down.clicked.connect(lambda: self.update_gratings('down'))
        self.pushButton_left.clicked.connect(lambda: self.update_gratings('left'))
        self.pushButton_right.clicked.connect(lambda: self.update_gratings('right'))

    def update_gratings(self, direction):
        step = float(self.lineEdit_step.text())
        grating_x = float(self.gratingx_lineEdit.text())
        grating_y = float(self.gratingy_lineEdit.text())
        if direction == 'center':
            self.gratingx_lineEdit.setText(str(0))
            self.gratingy_lineEdit.setText(str(0))
        elif direction == 'up':
            self.gratingy_lineEdit.setText(str(grating_y + step))
        elif direction == 'down':
            self.gratingy_lineEdit.setText(str(grating_y - step))
        elif direction == 'left':
            self.gratingx_lineEdit.setText(str(grating_x + step))
        elif direction == 'right':
            self.gratingx_lineEdit.setText(str(grating_x - step))

    def get_params(self):
        grating_x = float(self.gratingx_lineEdit.text())
        grating_y = float(self.gratingy_lineEdit.text())
        return grating_x, grating_y


class astigmatismUi(QtWidgets.QWidget):
    def __init__(self, slm_gui):
        super(astigmatismUi, self).__init__()
        uic.loadUi(os.path.join(os.path.dirname(__file__), 'ui_astigmatism.ui'), self)
        self.slm_gui = slm_gui
        self._connect()

    def _connect(self):
        self.pushButton_center.clicked.connect(lambda: self.update_astigmatism('center'))
        self.pushButton_up.clicked.connect(lambda: self.update_astigmatism('up'))
        self.pushButton_down.clicked.connect(lambda: self.update_astigmatism('down'))
        self.pushButton_left.clicked.connect(lambda: self.update_astigmatism('left'))
        self.pushButton_right.clicked.connect(lambda: self.update_astigmatism('right'))

    def update_astigmatism(self, direction):
        step = float(self.lineEdit_step.text())
        astigmatism_x = float(self.astigmatismx_lineEdit.text())
        astigmatism_y = float(self.astigmatismy_lineEdit.text())
        if direction == 'center':
            self.astigmatismx_lineEdit.setText(str(0))
            self.astigmatismy_lineEdit.setText(str(0))
        elif direction == 'up':
            self.astigmatismy_lineEdit.setText(str(astigmatism_y + step))
        elif direction == 'down':
            self.astigmatismy_lineEdit.setText(str(astigmatism_y - step))
        elif direction == 'left':
            self.astigmatismx_lineEdit.setText(str(astigmatism_x + step))
        elif direction == 'right':
            self.astigmatismx_lineEdit.setText(str(astigmatism_x - step))

    def get_params(self):
        astigmatism_x = float(self.astigmatismx_lineEdit.text())
        astigmatism_y = float(self.astigmatismy_lineEdit.text())
        return astigmatism_x, astigmatism_y


class focusUi(QtWidgets.QWidget):
    def __init__(self, slm_gui):
        super(focusUi, self).__init__()
        uic.loadUi(os.path.join(os.path.dirname(__file__), 'ui_focus.ui'), self)
        self.slm_gui = slm_gui
        self._connect()

    def _connect(self):
        # Connects the offset slider to the lineEdits
        self.lineEdit_step.returnPressed.connect(self.update_lineedit)
        self.lineEdit_offset.returnPressed.connect(self.update_lineedit)
        self.slider.valueChanged.connect(self.update_lineedit)
        self.lineEdit_value.returnPressed.connect(self.update_slider)
        # self.offset_slider.valueChanged.connect(self.make)

    def update_lineedit(self):
        step_size = float(self.lineEdit_step.text())
        offset = float(self.lineEdit_offset.text())
        steps = self.slider.value()
        value = offset + steps * step_size

        self.lineEdit_value.setText(str(value))

    def update_slider(self):
        value = float(self.lineEdit_value.text())
        step_size = float(self.lineEdit_step.text())
        offset = float(self.lineEdit_offset.text())

        steps = int((value - offset) / step_size)
        self.slider.setValue(steps)

    def get_params(self):
        curvature = float(self.lineEdit_value.text())
        return curvature,


class vortexbeamUi(QtWidgets.QWidget):
    def __init__(self, slm_gui):
        super(vortexbeamUi, self).__init__()
        uic.loadUi(os.path.join(os.path.dirname(__file__), 'ui_vortexbeam.ui'), self)
        self.slm_gui = slm_gui

    def get_params(self):
        order = int(float(self.lineEdit_order.text()))
        angle = float(self.lineEdit_angle.text())
        return order, angle


class linear_lutUi(QtWidgets.QWidget):
    def __init__(self, slm_gui):
        super(linear_lutUi, self).__init__()
        uic.loadUi(os.path.join(os.path.dirname(__file__), 'ui_linear_lut.ui'), self)
        self.slm_gui = slm_gui
        self._connect()

    def _connect(self):
        # Connects the offset slider to the lineEdits
        self.offset_lineEdit_step.returnPressed.connect(self.update_offset_lineedit)
        self.offset_lineEdit_offset.returnPressed.connect(self.update_offset_lineedit)
        self.offset_slider.valueChanged.connect(self.update_offset_lineedit)
        self.offset_lineEdit.returnPressed.connect(self.update_offset_slider)
        self.offset_slider.valueChanged.connect(self.slm_gui.make)

        # Connects the contrast slider to the lineEdits
        self.contrast_lineEdit_step.returnPressed.connect(self.update_contrast_lineedit)
        self.contrast_lineEdit_offset.returnPressed.connect(self.update_contrast_lineedit)
        self.contrast_slider.valueChanged.connect(self.update_contrast_lineedit)
        self.contrast_lineEdit.returnPressed.connect(self.update_contrast_slider)
        self.contrast_slider.valueChanged.connect(self.slm_gui.make)

    def update_offset_lineedit(self):
        step_size = float(self.offset_lineEdit_step.text())
        offset = float(self.offset_lineEdit_offset.text())
        steps = self.offset_slider.value()
        value = offset + steps * step_size

        self.offset_lineEdit.setText(str(value))

    def update_offset_slider(self):
        value = float(self.offset_lineEdit.text())
        step_size = float(self.offset_lineEdit_step.text())
        offset = float(self.offset_lineEdit_offset.text())

        steps = int((value - offset) / step_size)
        self.offset_slider.setValue(steps)

    def update_contrast_lineedit(self):
        step_size = float(self.contrast_lineEdit_step.text())
        offset = float(self.contrast_lineEdit_offset.text())
        steps = self.contrast_slider.value()
        value = offset + steps * step_size

        self.contrast_lineEdit.setText(str(value))

    def update_contrast_slider(self):
        value = float(self.contrast_lineEdit.text())
        step_size = float(self.contrast_lineEdit_step.text())
        offset = float(self.contrast_lineEdit_offset.text())

        steps = int((value - offset) / step_size)
        self.contrast_slider.setValue(steps)

    def get_params(self):
        contrast = float(self.contrast_lineEdit.text())
        offset = float(self.offset_lineEdit.text())
        return contrast, offset
