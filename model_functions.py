# -*- coding: utf-8 -*-
"""
Miscellaneous model functions
"""
import numpy as np
import parameters as p
from scipy import interpolate

#%%
''' Excitatory and inhibitory input functions P(t) and Q(t)'''
# Returns P, the external input given to excitatory cells, used in arg_e in algebraic_variables.py
# t: time vector
# input_data: experimental data used for the neuronal input functions
# location: for use inside a function ('in function') or an array for use outside of a function ('out function')
def input_P(t, input_data, location):  
    
    if p.InputCase == 'SquarePulse' or p.InputCase == 'Inhibitory': # Square pulse
        startpulse_end = p.startpulse + p.lengthpulse
        P = p.Pinput * (np.heaviside(t - p.startpulse,1) - np.heaviside(t - startpulse_end,1))
        
    elif p.InputCase == 'ZhengData' or p.InputCase == 'ZhengFittedParams': # Use neural input data        
        if location == 'in function':     
            if t/p.dt < p.Tend/p.dt:
                index_t = int(t/p.dt) # find corresponding index of the input_data array
            else:
                index_t = int(p.Tend/p.dt-1)
            P = p.Pinput * input_data[index_t]
        else: # for generating plotting vector
            P = p.Pinput * input_data  
            
    elif p.InputCase == '2SquarePulses': # Square pulse followed by second pulse
        startpulse_end = p.startpulse + p.lengthpulse;
        secondpulse_end = p.secondpulse + p.secondlength;
        P = p.Pinput * (np.heaviside(t - p.startpulse,1) - np.heaviside(t - startpulse_end,1) + np.heaviside(t - p.secondpulse,1) - np.heaviside(t - secondpulse_end,1));  
    return P

# Returns Q, the external input given to inhibitory cells, used in arg_i in algebraic_variables.py
# t: time vector
# input_data: experimental data used for the neuronal input functions
# location: for use inside a function ('in function') or an array for use outside of a function ('out function')
def input_Q(t, input_data, location):   
    
    if p.InputCase == 'SquarePulse' or p.InputCase == 'Inhibitory':  # Square pulse
        startpulse_end = p.startpulse + p.lengthpulse
        Q = p.Qinput * (np.heaviside(t - p.startpulse,1) - np.heaviside(t - startpulse_end,1))
        
    elif p.InputCase == 'ZhengData' or p.InputCase == 'ZhengFittedParams': # Use neural input data
        if location == 'in function':
            if t/p.dt < p.Tend/p.dt:
                index_t = int(t/p.dt) # find corresponding index of the input_data array
            else:
                index_t = int(p.Tend/p.dt-1)
            Q = p.Qinput * input_data[index_t]
        else: # for generating plotting vector
            Q = p.Qinput * input_data  
            
    elif p.InputCase == '2SquarePulses': # Square pulse followed by second pulse
        startpulse_end = p.startpulse + p.lengthpulse;
        secondpulse_end = p.secondpulse + p.secondlength;
        Q = p.Qinput * (np.heaviside(t - p.startpulse,1) - np.heaviside(t - startpulse_end,1) + np.heaviside(t - p.secondpulse,1) - np.heaviside(t - secondpulse_end,1));  
    return Q

''' Excitatory and inhibitory input function T(t)'''
# Returns T, the external input (thalamic) to both excitatory and inhibitory cells based on Pinto et al.
# t: time vector
# input_data: experimental data used for the neuronal input functions
# location: for use inside a function ('in function') or an array for use outside of a function ('out function')
# total_waveform: an array containing the fast triangular input pulses
def input_T(t, input_data, total_waveform, location): 
    
    if p.InputCase == 'ThalamicTrianglesZheng':       # Multiple fast triangular input scaled by Zheng data
        if location == 'in function':
            if t/p.dt < p.Tend/p.dt:
                index_t = int(t/p.dt) # find corresponding index based on time t
            else:
                index_t = int(p.Tend/p.dt-1)
            T = p.Tinput * total_waveform[index_t] * input_data[index_t]
        else: # for generating plotting vector outside of function
            T = p.Tinput * total_waveform * input_data  
            
    elif p.InputCase == 'ThalamicTriangles':     # Multiple fast triangular input pulses
        if location == 'in function':
            if t/p.dt < p.Tend/p.dt:
                index_t = int(t/p.dt) # find corresponding index based on time t
            else:
                index_t = int(p.Tend/p.dt-1)
            T = p.Tinput * total_waveform[index_t]
        else: # for generating plotting vector outside of function
            T = p.Tinput * total_waveform          
        
    elif p.InputCase == 'ThalamicSquarePulse':     # Simple input block
        startpulse_end = p.startpulse + p.lengthpulse
        T = p.Tinput * (np.heaviside(t - p.startpulse,1) - np.heaviside(t - startpulse_end,1))
        
    return T  

# Creates an array of the fast triangular pulses, returns total_waveform which is used in input_T function
def T_pulses():
    # Thalamic input: multiple fast triangular pulses (total_waveform) scaled by Zheng data (interp_neural_wh)
    if p.InputCase == 'ThalamicTriangles' or p.InputCase == 'ThalamicTrianglesZheng': # Model with InputCase (Pinto et al) input - triangles for thalamic input T(t), overwrites input_data
        # Parameters
        Hz = 5                   # [Hz or s**-1]                              # [individualpulses s**-1]
        peak_time = 1            # [ms]                                      # Thalamic input increase time
        
        # Expressions
        msperpulse = 1000 / Hz                                            # [ms individualpulse**-1]
        index_peak_time = peak_time + 1                                    # Individual pulse increase length
        down_time = msperpulse - index_peak_time #msperpulse - index_peak_time                           # Individual pulse decrease length [ms], if we want the triangles to last only 15ms then use down_time = 15 - index_peak_time
        rest_time = msperpulse - (down_time + index_peak_time)             # Rest time individual pulse
        numberOfPeriods = int(p.lengthpulse / msperpulse)                      # Number of individual pulses in neuronal stim

        # Make start and end time resting values
        start_time = p.startpulse
        end_time = p.Tend - (p.startpulse + p.lengthpulse) 
        base_input = 0.0  # Base thalamic input
        startSignal = np.linspace(base_input, base_input, int(start_time/p.dt))
        endSignal = np.linspace(base_input, base_input, int(end_time/p.dt))

        # Construct one individual pulse cycle     
        firstamp = 1.07  # Thalamic input amplitude
        firstRiseSignal = np.linspace(base_input, firstamp, int(index_peak_time))
        firstFallingSignal = np.linspace(firstamp, base_input, int(down_time))
        restSignal = np.linspace(base_input, base_input, int(rest_time))
        firstCycle = np.append(firstRiseSignal, firstFallingSignal[1:-1])
        firstCycle = np.append(firstCycle, restSignal)

        # Now replicate this cycle several (numberOfPeriods) times.
        waveform = np.tile(firstCycle, (1, numberOfPeriods)).flatten()

        # Interpolate so that it matches with the timestep p.dt        
        waveform_interp_f = interpolate.interp1d(np.linspace(0, p.lengthpulse, len(waveform)), waveform)
        waveform_interp = waveform_interp_f(np.arange(0, p.lengthpulse+p.dt, p.dt))  
        
        # Add the bits before and after the stimulus
        total_waveform = np.append(startSignal, waveform_interp)
        total_waveform = np.append(total_waveform, endSignal)
    else:
        total_waveform = []
    
    return total_waveform
