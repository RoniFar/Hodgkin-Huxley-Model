import numpy as np

class ModelHH:
    """the model consist all the paramaters and calaulation methods
    parameter:
    * EK=-12    Potassium potential (mV)
    * ENa=115   Sodium potential (mV)
    * El=10.6   Leak potential (mV)
    * gK=36     Average potassium channel conductance per unit area (mS/cm^2)
    * gNa=120   Average sodium channel conductance per unit area (mS/cm^2)
    * gl=0.3    Average leak channel conductance per unit area (mS/cm^2)
    * Cm=1      Membrane capacitance per unit area (uF/cm^2). between 0.7 and 1.
    * Vm=0      Membrane potential (mV)
    * IK=0      Potassium current (uA/cm^2)
    * INa=0     Sodium current (uA/cm^2)
    * Il=0      Leak current (uA/cm^2)
    """

    class Gate:
        "the kineticts of the channels and the state (open)"
        def __init__(self, alpha=0, beta=0, state=0, tau=0, inf=0):
            self.alpha = alpha
            self.beta = beta
            self.state = state
            self.tau, self.inf = tau, inf

        def AlphaBetaUpdate(self, Vm, gate):
            """update channel's kinetics and gate time constant according to formula
            in steady state"""
            if gate == "n":
                self.alpha = (0.01 * (10 - Vm)) / (np.exp((10 - Vm) / 10) - 1)
                self.beta = 0.125 * np.exp(-Vm / 80)
            if gate == "m":
                self.alpha = (0.1 * (25 - Vm)) / (np.exp((25 - Vm) / 10) - 1)
                self.beta = 4 * np.exp((-Vm) / 18)
            elif gate == "h":
                self.alpha = 0.07 * np.exp(-Vm / 20)
                self.beta = 1 / (np.exp((30 - Vm) / 10) + 1)
            self.setTauInf()

        def setTauInf(self):
            self.tau = 1 / (self.alpha + self.beta)
            self.inf = self.alpha / (self.alpha + self.beta)

        def dGdt(self, Vm, gate):
            "calculate by the equ dG/dt = (G_inf(V)-G)/tauG(V), where G in {n,m,h}"
            self.AlphaBetaUpdate(Vm, gate)
            return (self.inf - self.state)/self.tau

    def __init__(self, EK=-12, ENa=115, El=10.6, gK=36, gNa=120, gl=0.3, Cm=1, Vm=0, \
                 IK=0, INa=0, Il=0, s_a=1):
        """define all model pamaters and components"""
        self.EK, self.ENa, self.El = EK, ENa, El
        self.gK, self.gNa, self.gl = gK * s_a, gNa * s_a, gl * s_a
        self.n, self.m, self.h = self.Gate(), self.Gate(), self.Gate()
        self.Cm = Cm    # membrane cappacity
        self.Vm = Vm    # Resting state membrane potential , starting voltag
        self.IK, self.INa, self.Il = IK, INa, Il

    def NumericUpdate(self, dt, Iinj):
        """ Calculating by Euler formula to numeric process:
            dy/dt = f(y,t)
            y(t+∆t) ~= y(t) + (dy(t)/dt) * ∆t
            the function dGdt calculate the value of dy(t)/dt"""
        self.n.state += self.n.dGdt(self.Vm, "n") * dt
        self.m.state += self.m.dGdt(self.Vm, "m") * dt
        self.h.state += self.h.dGdt(self.Vm, "h") * dt
        self.Vm += ((self.dVdt() + Iinj) * dt) / self.Cm

    def dVdt(self):
        """calculate ion current and update"""
        self.IK = self.gK * np.power(self.n.state, 4) * (self.EK - self.Vm)
        self.INa = self.gNa * np.power(self.m.state, 3) * self.h.state * (self.ENa -self.Vm)
        self.Il = self.gl * (self.El - self.Vm)
        return self.IK + self.INa + self.Il


class HHcompartment:
    """Hodgkin-Huxley compartments model
        * Cm=1      Membrane capacitance per unit area (uF/cm^2). between 0.7 and 1.
        * ro=r_a*a/L   Specific resistance of axoplasm (ohm*cm)
        * a=0.00005    Radius of axon (cm)
        * L=0.000010    Length of axon (cm)
        * u=0.1880   Velocity of moving axes (m/sec)
        * a=0.0239     Radius (cm)
        * r_a=35.4  Exoplasmic resistivity (ohm*cm)
        * g_a = 1/r_a
        * g_ij=0    resistance coupling between neighboring compartmants from i to j
        note: where the unit are µm we need to mo"""

    def __init__(self, Num=2, L=0.055, a=0.0239, ro=35.4, Cm=1):
        self.Num = Num
        self.L = L
        self.a = a
        self.ro = ro
        self.surfaceArea = 2 * np.pi * self.a * self.L
        self.crossSection = np.pi * np.power(self.a, 2)
        self.r_a = self.ro * (self.L / self.crossSection)
        self.g_a = 1 / self.r_a
        self.C = Cm * self.surfaceArea

    def currBetween(self, idx):
        """didn't use but the logic here is numeric for PDE from second order"""
        if idx == 0:
            I_i = self.g_a * (self.comp[1].Vm - self.comp[0].Vm)
        elif idx == self.Num - 1:
            I_i = self.g_a * (self.comp[idx-1].Vm - self.comp[idx].Vm)
        else:
            I_i = self.g_a * (self.comp[idx+1].Vm + self.comp[idx-1].Vm - 2 * self.comp[idx].Vm)
        return I_i








