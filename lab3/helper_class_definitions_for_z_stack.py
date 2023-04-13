class sub_SMC100_move(QRunnable):
    def __init__(self, SMC100, axis, value):
        """ axis is a string and value is absolute position """
        super().__init__()
        self.SMC100 = SMC100
        self.axis = axis
        self.value = value        
        
    def run(self):
        global stage_in_use
        if stage_in_use ==True:
            print('Wait for stage release')
        while stage_in_use ==True:            
            time.sleep(0.5)
        print('Moving SMC100')
        stage_in_use = True
        pos_stage = self.SMC100.get_position(self.axis)
        goto_pos = float(pos_stage[0]) + self.value
        self.SMC100.move(goto_pos, self.axis, relative=False)
        pos_stage = self.SMC100.get_position(self.axis)
        print('Position axis ' + str(self.axis) + ': ' + str(pos_stage) + ' mm')
        stage_in_use = False 
        