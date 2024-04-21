import numpy as np
from functions import *
from HHtypes import *

# for saving the figs change the val of "savefig" to True
savefig, showfig = True, False

# time parameters
delta = 1e-03
maxTime = 100
tVec = np.arange(0, maxTime+delta, delta)           # adding delta to maxTime to make sure the vec ends with maxTime val

# Q1
HH1 = ModelHH()                                     # HH model
amp = 10                                            # Current amplitude in µA
StartInj, EndInj = 5, 10                            # Current injection for short time
# StartInj, EndInj = 5, 95                          # Current injection for long time
Iinj = Iinject(tVec, StartInj, EndInj, amp)     # creating vector of inject current
Data1 = sim(tVec)                                   # experiment model
Data1.IVupdate(HH1, delta, Iinj)
Data1.APPT()

# Q1 plot
plot_q1(amp, StartInj, EndInj, Iinj, tVec, maxTime, Data1, savefig=savefig, show=showfig, fignum=1)

# Q2
# creating the model
compNum = 10
AxonParam = HHcompartment()
HH_lst = np.empty((0,compNum))
Data_lst = np.empty((0,compNum))
HH_lst = np.append(HH_lst, HH1)
Data_lst = np.append(Data_lst, Data1)

# iteration over compartments
HH_lst, Data_lst = CompCalc(compNum, AxonParam, tVec, delta, HH_lst, Data_lst)

# Q2 plot of action potential progress
plot_q2_ap_prog(amp, StartInj, EndInj, Iinj, tVec, maxTime, Data_lst, compNum, savefig=savefig, show=showfig, fignum=1)

# changing the radius
# calculate the speed between neighboring compartments to make the code run faster.
compNum = 3
measureNum = 15         # make smaller for fast run
a_min = AxonParam.a / 2
a_max = AxonParam.a * (3 / 2)
delta_a = (a_max - a_min) / measureNum
a_list = np.arange(a_min, a_max, delta_a)
ConductSpeed = np.empty(measureNum)

HH_iter_lst = np.empty((0,compNum))
Data_iter_lst = np.empty((0,compNum))
HH_iter_lst = np.append(HH_iter_lst, HH1)
Data_iter_lst = np.append(Data_iter_lst, Data1)

for ii in range(measureNum):
    ConductSpeed[ii] = calc_speed(a_list[ii], compNum, tVec, delta, HH_iter_lst, Data_iter_lst)

plot_q2_radius_speed(a_list, ConductSpeed, savefig=savefig, show=showfig, fignum=1)




