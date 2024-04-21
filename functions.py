import numpy as np
import os
import matplotlib.pyplot as plt
from HHtypes import *
from Experiments import Current as sim


def Iinject(t, injStart, injEnd, amp):
    """ this function creating an array I, that represents the current step """
    I = np.zeros(np.shape(t))
    S = t >= injStart
    E = t <= injEnd
    mask = S == E   # creating time window for injection in mili-sec
    I[mask] = amp   # amplitude [micro-amper]
    return I


def CompCalc(compNum, AxonParam, tVec, delta, hh_lst, data_lst, only_picks=False):
    for i in range(1, compNum):
        hh_lst = np.append(hh_lst, ModelHH(s_a=AxonParam.surfaceArea, Cm=AxonParam.C))
        Data = sim(tVec)
        Data.compUpdate(hh_lst[i], delta, data_lst[i - 1].VmVec, AxonParam)
        Data.APPT()
        data_lst = np.append(data_lst, Data)
    if only_picks:
        picks = np.empty(compNum)
        for ii in range(compNum):
            picks[ii] = data_lst[ii].APpicks
        return picks
    return hh_lst, data_lst


def calc_speed(A, compNum, tVec, delta, HH_iter_lst, Data_iter_lst):
    AxonParam_iter = HHcompartment(a=A)
    Picks = CompCalc(compNum, AxonParam_iter, tVec, delta, HH_iter_lst, Data_iter_lst, only_picks=True)
    # distance over time equal velocity
    return AxonParam_iter.L / (Picks[-1] - Picks[1])

########################################################################################################################
################################################ Ploting ###############################################################
########################################################################################################################
def plot_q1(amp, StartInj, EndInj, Iinj, tVec, maxTime, Data1, savefig=False, show=False, fignum=1):
    """plot function for Q1.
    savefig geta True/False"""
    dirpath = os.path.dirname(os.getcwd())
    fold_path = dirpath + '/results/'

    xjump = 10
    fig, (ax1, ax2) = plt.subplots(2, 1, sharex=True, figsize=(12, 10))
    plt.suptitle('Hodgkin-Huxley Model \n Steady Current injection of ' \
                 + str(amp) + '[µA] for ' + str(EndInj - StartInj) + '[mS]', fontsize=15)
    ax1.plot(tVec, Data1.VmVec)
    ax1.set_title('Voltage for time curve', fontsize=10)
    ax1.set_ylim(-20, 120)
    ax1.set_ylabel("Membrane potential [mV]")

    ax2.plot(tVec, Iinj, 'y')
    ax2.set_title('Current injection for time curve', fontsize=10)
    ax2.plot(tVec, np.zeros(len(tVec)), color='black', linestyle='dashed', linewidth=0.5)
    ax2.set_ylim(-5, 30)
    ax2.set_ylabel("Injection Current [µA]")
    ax2.set_xlabel("time [mS]")
    ax2.set_xticks(np.arange(0, maxTime + xjump, xjump))
    if savefig:
        if not os.path.isdir(fold_path):
            os.makedirs(fold_path)
        figname = "CurrentInjection_" + str(amp) + "_duration_" + str(EndInj - StartInj) + '_' + str(fignum)
        plt.savefig(fold_path + figname)
        print('figured - "' + figname + '" was saved in ' + fold_path)
    if show:
        plt.show()

def plot_q2_ap_prog(amp, StartInj, EndInj, Iinj, tVec, maxTime, Data_lst, compNum, savefig=False, show=False, fignum=1):
    dirpath = os.path.dirname(os.getcwd())
    fold_path = dirpath + '/results/'

    xjump = 10
    fig, (ax1, ax2) = plt.subplots(2, 1, sharex=True, figsize=(12, 10))
    plt.suptitle('Hodgkin-Huxley Compartments Model \n Steady Current injection of ' \
                 + str(amp) + '[µA] for ' + str(EndInj - StartInj) + '[mS]')
    for i in range(compNum):
        ax1.plot(tVec, Data_lst[i].VmVec, label="compartment " + str(i+1))
    ax1.set_title('Voltage for time curve', fontsize=10)
    ax1.set_ylabel("Membrane potential [mV]")
    ax1.legend()

    ax2.plot(tVec, Iinj, 'y')
    ax2.set_title('Current injection for time curve', fontsize=10)
    ax2.plot(tVec, np.zeros(len(tVec)), color='black', linestyle='dashed', linewidth=0.5)
    ax2.set_ylim(-5, 30)
    ax2.set_ylabel("Injection Current [µA]")
    ax2.set_xlabel("time [mS]")
    ax2.set_xticks(np.arange(0, maxTime + xjump, xjump))
    if savefig:
        if not os.path.isdir(fold_path):
            os.makedirs(fold_path)
        figname = 'CompartmentsModelProgress_' + str(fignum)
        plt.savefig(fold_path + figname)
        print('figured - "' + figname + '" was saved in ' + fold_path)
    if show:
        plt.show()


def plot_q2_radius_speed(radius_lst, conduct_speed, savefig=False, show=False, fignum=1):
    dirpath = os.path.dirname(os.getcwd())
    fold_path = dirpath + '/results/'

    plt.subplots(figsize=(12, 7))
    plt.suptitle('conductance speed')
    plt.plot(radius_lst, conduct_speed)
    plt.xlabel("radius [cm]")
    plt.ylabel("conductance speed[mS]")
    if savefig:
        if not os.path.isdir(fold_path):
            os.makedirs(fold_path)
        figname = 'CompartmentsModel_Radius_vs_Speed_' + str(fignum)
        plt.savefig(fold_path + figname)
        print('figured - "' + figname + '" was saved in ' + fold_path)
    if show:
        plt.show()