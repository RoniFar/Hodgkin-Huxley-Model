import numpy as np

class Current:

    def __init__(self, time):
        self.tVec = time
        self.VmVec = np.zeros(len(time))
        self.CurrVec = np.zeros(len(time))
        self.APpicks = 0

    def IVupdate(self, model, dt, Iinj):
        """:parameter
        model = class that has the functions and parameters of hodgkin-huxley model
        dt = int time interval between update
        Iinj = ndarray with len of tVec that contain the current amplitude at time t"""
        model.NumericUpdate(dt, Iinj[0])  # initialized the model parameters for steady state
        self.VmVec[0] = model.Vm
        model.n.state, model.m.state, model.h.state = model.n.inf, model.m.inf, model.h.inf     # for loop update
        for t in range(1, len(self.tVec)):
            model.NumericUpdate(dt, Iinj[t])
            self.VmVec[t] = model.Vm

    def compUpdate(self, model, dt, V_prev_comp, axonparam):
        model.NumericUpdate(dt, V_prev_comp[0])  # initialized the model parameters for steady state
        self.VmVec[0] = model.Vm
        model.n.state, model.m.state, model.h.state = model.n.inf, model.m.inf, model.h.inf     # for loop update
        self.CurrVec[0] = 0

        for t in range(1, len(self.tVec)):
            self.CurrVec[t] = (V_prev_comp[t] - self.VmVec[t-1]) / axonparam.r_a
            model.NumericUpdate(dt, self.CurrVec[t])
            self.VmVec[t] = model.Vm

    def APPT(self):
        """action potential picks time"""
        idx = self.VmVec.argmax()
        self.APpicks = self.tVec[idx]



