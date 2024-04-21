# Hodgkin-Huxley-Model
This code calaulate cell voltage using Hodgkin–Huxley equations for giant squid axon.

For running the code please run "main.py".
    NOTE: If you want to save the visualization change "savefig" variable in line 6 to True
    time parameters: define "delta" as the time interval between 2 measurements.
                     define "maxTime" as the measurement full interval.
    NOTE: if you want the code to run faster change "measureNum" in line 44 to smaller number.

    The main script will call all the other necessary scrips for calculation and visualization


"HHtypes.py" oncluded all the classes of the model:
    # ModelHH - creating the Hodgkin-Huxley Model and includes functions for the numeric updates
    # HHcompartment - class for axons parameters, default for Giant squid axon.

"Experiments.py" creating class with all the necessary functions to calculate the current in the system, and saving the
    data in each time measure about the current, voltage and the first pick of the neuron (the pick of the first
    action potential).

"function.py" organizes all the scrips that in multiple uses as functions, and all the printing command that
    does not necessary for the logic of model and code.
