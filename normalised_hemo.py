# -*- coding: utf-8 -*-
"""
Normalised BOLD and hemodynamics
"""

import numpy as np
import parameters as p

# Returns the class 'a' with added normalised variables that depend on the steady state conditions
# a: class containing the algebraic variables
# v: class containing the state variables
# t: time vector
def solve_normalised_hemodynamics(a,v,t):
    
    # Find the index for a time right before the start of stimulation
    if p.startpulse < p.Tend:
        preNeuronalStimTime = np.where(t == p.startpulse) 
        preNeuronalStimTime = int(preNeuronalStimTime[0][0])
    else: # 1 second before the end of the simulation (if there is no stimulation within the simulation period)
        preNeuronalStimTime = np.where(t == p.Tend - 1e3) 
        preNeuronalStimTime = int(preNeuronalStimTime[0][0])
    
    # Find steady state conditions (just before the start of stimulation)
    CBF_0 = a.CBF[preNeuronalStimTime]
    CBV_0 = v.CBV[preNeuronalStimTime]
    HBR_0 = v.HbR[preNeuronalStimTime]
    CMRO2_0 = a.CMRO2[preNeuronalStimTime]
    
    # Normalised variables
    CBF_N = a.CBF/CBF_0         # Cerebral blood flow
    CBV_N = v.CBV/CBV_0         # Cerebral blood volume
    HBR_N = v.HbR/HBR_0         # Deoxyhemoglobin concentration
    CMRO2_N = a.CMRO2/CMRO2_0   # CMRO2
    
    # Solve the other variables
    HBT_N = CBF_N * HBR_N / CMRO2_N                                         # Total hemoglobin (normalised), equivalent to HBR_N .* p.E_0 ./ OEF
    HBO_N = (HBT_N - 1) - (HBR_N - 1) + 1                                   # Oxyhemoglobin (normalised)
    BOLD_N = 100 * p.V_0 * ( p.a_1 * (1 - HBR_N) - p.a_2 * (1 - CBV_N) )    # BOLD (percentage increase from 0)
    
    # Add these to the class of algebraic variables
    a.CBF_N = CBF_N
    a.CBV_N = CBV_N
    a.HBR_N = HBR_N
    a.CMRO2_N = CMRO2_N
    a.HBT_N = HBT_N
    a.HBO_N = HBO_N
    a.BOLD_N = BOLD_N
    
    return a