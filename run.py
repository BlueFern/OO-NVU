# -*- coding: utf-8 -*-
"""
NVU 3.2?

Run the NVU model code from here using ODEINT to solve the differential equations
Also plot a selection of variables
"""

'''
Notes (for Spyder IDE): 
    - to have plots open in new windows rather than in the console, go to Preferences->IPython Console->Graphics and set Graphics Backend to "Automatic"
    - to show all variables in Variable explorer (including figures and classes), go to Preferences->Variable Explorer and turn off "Exclude unsupported data types" 
    - to allow for plotting over existing figures, go to Preferences->Run and turn off "Remove all variables before execution"
    - restart Spyder to take these into effect
'''

# Import modules
from scipy.integrate import odeint
import numpy as np
import time
import matplotlib.pyplot as plt
import warnings
from matplotlib.patches import Rectangle

# Local files
from indices import set_indices
from ICs import set_initial_conditions
from ODEsystem import func
from algebraic_variables import set_algebraic_variables
from state_variable_names import initialise_var_names
from normalised_hemo import solve_normalised_hemodynamics
from plotting_functions import plot_variables_singles, plot_variables_samegraph
from model_functions import T_pulses
import parameters as p
import import_mat_files as im

# Suppress runtime warnings, comment this out if something is going wrong
warnings.filterwarnings("ignore", category=RuntimeWarning)

# Calculate the start time
start_time = time.time()

# Initialise indices, initial conditions and time vector
idx = set_indices()
u0 = set_initial_conditions(idx)
t = np.arange(start=0, stop=p.Tend+p.dt, step=p.dt)

# Import the experimental neural input profile from mat file, for InputCase = 'ZhengData', 'ThalamicTrianglesZheng' or 'ZhengFittedParams'
input_data = im.import_Zheng_neural_data(t)

# Generate the multiple triangular pulse input profile for the thalamic input T(t), for InputCase = 'ThalamicTriangles' or 'ThalamicTrianglesZheng'
total_waveform = T_pulses()

# Solve ODE system
print("\n Model simulation successfully started\n")
u, solver_details = odeint(func, u0, t, args=(idx,input_data,total_waveform,start_time), hmax=1e2, full_output=1)
if solver_details['message'] == 'Integration successful.':
    del solver_details

# Print total time taken
solve_time = time.time() - start_time
if solve_time < 60.0:
    print("\n\n--- Took {0:.0f} seconds to solve ---".format(solve_time))
else:
    minutes, seconds= divmod(solve_time, 60)
    print("\n\n--- Took {0:.0f} minutes and {1:.0f} seconds to solve ---".format(minutes, seconds))
del start_time, solve_time
    
# Extract the state variables for easy plotting etc
v = initialise_var_names(u, idx, 'out function')

# Extract the algebraic variables for easy plotting etc
a = set_algebraic_variables(u, t, idx, input_data, total_waveform, 'out function')
del total_waveform, input_data

# Solve for the normalised CBF, CBV, CMRO2, HbO, HbR, HbT and BOLD and put into 'a' (solved separately from odeint as they depend on the steady state conditions)
a = solve_normalised_hemodynamics(a,v,t)

# The variables that can be plotted sorted into alphabetical order
state_vars = sorted(vars(v).keys())
#print("\nState variables: \n", state_vars, "\n")
alg_vars = sorted(vars(a).keys())
#print("Algebraic variables: \n", alg_vars,"\n")

'''Plotting parameters'''
if p.startpulse < p.Tend:
    time = t/1e3 - p.startpulse/1e3     # Time in seconds, normalised so t=0 at the start of stimulation
    xlim1 = p.startpulse/1e3 - 10       # Set left xlimit to 10 sec before stimulation begins
else:
    time = t/1e3                        # Time in seconds, not normalised if there is no stimulation
    xlim1 = 0                           # Set left xlimit to 0 if there is no stimulation
xlim2 = p.Tend/1e3  # Set right xlimit to end of simulation
del t

# Set custom xlims if wanted
xlim1 = -1
xlim2 = 20

# %%
'''Plots using the automated plotting functions (in plotting_functions module) - just specify a list of things you want to plot from state and/or algebraic variables'''

## Plot all state variables in subplots
#fig_statevars = plot_variables_singles(time, a, v, xlim1, xlim2, plot_list=state_vars, title='All state variables', num=1)
#fig_statevars.savefig("./figures/State_variables.pdf", bbox_inches='tight')

## Plot all algebraic variables in subplots
#fig_algvars = plot_variables_singles(time, a, v, xlim1, xlim2, plot_list=alg_vars, title='All algebraic variables')

## Plot a single variable
#fig_singlevar = plot_variables_singles(time, a, v, xlim1, xlim2, plot_list=['O2'], title='Oxygen [mM]', figsize=(9,6), dpi=100)

### Plot a custom selection in subplots
plot_list = ['H_i','w_i','J_K_i','v_i','J_VOCC_i','Ca_i','R']    
fig_custom_subplots = plot_variables_singles(time, a, v, xlim1, xlim2, plot_list, alphabetical=0, num=14, show_stimulation=0)

## Plot a selection on the same graph
#plot_list = ['HBO_N','HBR_N','HBT_N']
#fig_custom_samegraph = plot_variables_samegraph(time, a, v, xlim1, xlim2, plot_list, title='Hemodynamics', legend=['HbO','HbR','HbT'], colourstyle=['r-','b-','g-'], figsize=(6,4), show_stimulation=1)
#plt.close(fig_custom_samegraph)

# %% 
'''Import and plot experimental data from file''' 

## CBF data from Zheng 2010
#fig_Zheng_CBF, cbf_time, cbf_Zheng = im.import_Zheng_cbf_data()

## Air and oxygen data from Berwick et al
# Choose whether to look at a specific region or all regions - 'All', 'Whisker', 'Vein', 'Artery', 'Parenchyma'. If left blank then 'All' is assumed
#fig_airOxyComparison, fig_airOxyTTest_HbO, fig_airOxyTTest_HbR, fig_airOxyTTest_HbT, airoxydata = im.import_Berwick_AirOxyData()
#plt.close(fig_airOxyComparison)
#plt.close(fig_airOxyTTest_HbO)
#plt.close(fig_airOxyTTest_HbR)
#plt.close(fig_airOxyTTest_HbT)

## LNAME data from Berwick et al
# Choose whether to look at a specific region or all regions - 'All', 'Whisker', 'Vein', 'Artery', 'Parenchyma'. If left blank then 'All' is assumed
#fig_LNAME_experiment, LNAME_data = im.import_Berwick_LNAME_7NI_data() # Import all data and save in a class
#time_LNAME, exp_HbO, exp_HbR, exp_HbT, exp_HbO_std, exp_HbR_std, exp_HbT_std = im.set_corresponding_LNAME_7NI_experimental_data(LNAME_data) # Find the experiment that corresponds to the current model simulation
#plt.close(fig_LNAME_experiment)

# %%
''' Custom section - do what you want here '''
'''Handy matplotlib guide: https://towardsdatascience.com/all-your-matplotlib-questions-answered-420dd95cb4ff'''
'''Use num = <number> in plt.figure() to plot on a specific figure - useful for plotting on existing figures'''

### Plot the air data, oxygen data, whisker normal data, and model results for HbO, HbR, HbT
#
#fig_airOxyHemoComparison = plt.figure(figsize=(13,5), dpi=80, edgecolor='k') 
##fig_airOxyHemoComparison.suptitle("Hshift = {0:.0f}".format(p.Hshift))
#
#fig_airOxyHemoComparison.add_subplot(1,3,1)
#plt.gca().get_yaxis().get_major_formatter().set_useOffset(False) # stop weird exponential notation
#plt.plot(time, a.HBO_N, label = 'HbO model')
#plt.plot(airoxydata.time_airOxy, airoxydata.HbOair_mean, ':', label = 'HbO air')
#plt.plot(airoxydata.time_airOxy, airoxydata.HbOoxy_mean, ':', label = 'HbO oxy')
##plt.plot(time_HbO, exp_HbO, '--', label = 'HbO exp')
#plt.title('HbO')
#plt.xlabel('Time [s]') 
#plt.xlim(-2,20)  
#plt.ylim(0.98, 1.09)  
#plt.legend(loc='best')
##plt.tight_layout()
#
#fig_airOxyHemoComparison.add_subplot(1,3,2)
#plt.gca().get_yaxis().get_major_formatter().set_useOffset(False) # stop weird exponential notation
#plt.plot(time, a.HBR_N, label = 'HbR model')
#plt.plot(airoxydata.time_airOxy, airoxydata.HbRair_mean, ':', label = 'HbR air')
#plt.plot(airoxydata.time_airOxy, airoxydata.HbRoxy_mean, ':', label = 'HbR oxy')
##plt.plot(time_HbR, exp_HbR, '--', label = 'HbR exp')
#plt.title('HbR')
#plt.xlabel('Time [s]') 
#plt.xlim(-2,20)  
#plt.ylim(0.91, 1.02)  
#plt.legend(loc='best')
##plt.tight_layout()
#
#fig_airOxyHemoComparison.add_subplot(1,3,3)
#plt.gca().get_yaxis().get_major_formatter().set_useOffset(False) # stop weird exponential notation
#plt.plot(time, a.HBT_N, label = 'HbT model')
#plt.plot(airoxydata.time_airOxy, airoxydata.HbTair_mean, ':', label = 'HbT air')
#plt.plot(airoxydata.time_airOxy, airoxydata.HbToxy_mean,  ':', label = 'HbT oxy')
##plt.plot(time_HbT, exp_HbT, '--', label = 'HbT exp')
#plt.title('HbT')
#plt.xlabel('Time [s]') 
#plt.xlim(-2,20)  
#plt.ylim(0.99, 1.05)  
#plt.legend(loc='best')
#plt.tight_layout()
#

#%%
#fig_airOxyStd = plt.figure(figsize=(14,5), dpi=80, edgecolor='k') 
#
#fig_airOxyStd.add_subplot(1,3,1)
#plt.gca().get_yaxis().get_major_formatter().set_useOffset(False) # stop weird exponential notation
#plt.errorbar(airoxydata.time_airOxy, airoxydata.HbOair_mean, airoxydata.HbOair_std, capsize=3, errorevery=12, elinewidth=1.5, capthick=1, label = 'Air')
#plt.errorbar(airoxydata.time_airOxy, airoxydata.HbOoxy_mean, airoxydata.HbOoxy_std, capsize=3, errorevery=12, elinewidth=1.5, capthick=1, label = 'Oxy')
#plt.title('HbO')
#plt.xlabel('Time [s]') 
#plt.xlim(-2,20)  
##plt.ylim(0.98, 1.09)  
#plt.legend(loc='best')
##plt.tight_layout()
#
#fig_airOxyStd.add_subplot(1,3,2)
#plt.gca().get_yaxis().get_major_formatter().set_useOffset(False) # stop weird exponential notation
#plt.errorbar(airoxydata.time_airOxy, airoxydata.HbRair_mean, airoxydata.HbRair_std, capsize=3, errorevery=12, elinewidth=1.5, capthick=1, label = 'Air')
#plt.errorbar(airoxydata.time_airOxy, airoxydata.HbRoxy_mean, airoxydata.HbRoxy_std, capsize=3, errorevery=12, elinewidth=1.5, capthick=1, label = 'Oxy')
#plt.title('HbR')
#plt.xlabel('Time [s]') 
#plt.xlim(-2,20)  
##plt.ylim(0.91, 1.02)  
#plt.legend(loc='best')
##plt.tight_layout()
#
#fig_airOxyStd.add_subplot(1,3,3)
#plt.gca().get_yaxis().get_major_formatter().set_useOffset(False) # stop weird exponential notation
#plt.errorbar(airoxydata.time_airOxy, airoxydata.HbTair_mean, airoxydata.HbTair_std, capsize=3, errorevery=12, elinewidth=1.5, capthick=1, label = 'Air')
#plt.errorbar(airoxydata.time_airOxy, airoxydata.HbToxy_mean, airoxydata.HbToxy_std, capsize=3, errorevery=12, elinewidth=1.5, capthick=1, label = 'Oxy')
#plt.title('HbT')
#plt.xlabel('Time [s]') 
#plt.xlim(-2,20)  
##plt.ylim(0.99, 1.05)  
#plt.legend(loc='best')

#%%

# Plot the hemodynamics against the whisker/opto experiment
#fig_hemocomparison = plt.figure(figsize=(6,4), dpi=80, edgecolor='k', num=3) 
#plt.gca().get_yaxis().get_major_formatter().set_useOffset(False) # stop weird exponential notation
#plt.plot(time, a.HBO_N, 'r', label = 'HbO')
#plt.errorbar(time_LNAME, exp_HbO, exp_HbO_std, capsize=3, errorevery=15, elinewidth=1, capthick=1, color='r', ls=':', label = 'HbO exp')
#plt.plot(time, a.HBR_N, 'b', label = 'HbR')
#plt.errorbar(time_LNAME, exp_HbR, exp_HbR_std, capsize=3, errorevery=15, elinewidth=1, capthick=1, color='b', ls=':', label = 'HbR exp')
#plt.plot(time, a.HBT_N, 'g', label = 'HbT')
#plt.errorbar(time_LNAME, exp_HbT, exp_HbT_std, capsize=3, errorevery=15, elinewidth=1, capthick=1, color='g', ls=':', label = 'HbT exp')
#plt.title('{0:s}, {1:s}'.format( p.NeuronType.capitalize(), p.NOswitch.capitalize()))
#plt.xlabel('Time [s]') 
#plt.xlim(-2,20)  
#plt.ylim(0.86, 1.08)  
#plt.legend(loc='best')
#plt.tight_layout()
##

# %%

### Plot a selection on subplots

#fig12 = plt.figure(figsize = (20,15), dpi=80, edgecolor='k', num = 2) 
#
#fig12.add_subplot(2,4,1)
#plt.gca().get_yaxis().get_major_formatter().set_useOffset(False) # stop weird exponential notation
#plt.plot(time, a.f_in)
#plt.title('f_{in}')
#plt.xlabel('Time [s]')
#plt.xlim(xlim1, xlim2)
#plt.ylim(0.99, 1.09)
##plt.tight_layout()
#
#fig12.add_subplot(2,4,2)
#plt.gca().get_yaxis().get_major_formatter().set_useOffset(False) # stop weird exponential notation
#plt.plot(time, a.OEF/p.E_0)
#plt.title('OEF/E_0')
#plt.xlabel('Time [s]')
#plt.xlim(xlim1, xlim2)
#plt.ylim(0.94, 1.01)
##plt.tight_layout()
#
#fig12.add_subplot(2,4,3)
#plt.gca().get_yaxis().get_major_formatter().set_useOffset(False) # stop weird exponential notation
#plt.plot(time, v.HbR)
#plt.title('HbR')
#plt.xlabel('Time [s]')
#plt.xlim(xlim1, xlim2)
#plt.ylim(0.96, 1.01)
##plt.tight_layout()
#
#fig12.add_subplot(2,4,4)
#plt.gca().get_yaxis().get_major_formatter().set_useOffset(False) # stop weird exponential notation
#plt.plot(time, v.CBV)
#plt.title('CBV')
#plt.xlabel('Time [s]')
#plt.xlim(xlim1, xlim2)
#plt.ylim(0.99, 1.03)
##plt.tight_layout()
#
#fig12.add_subplot(2,4,5)
#plt.gca().get_yaxis().get_major_formatter().set_useOffset(False) # stop weird exponential notation
#plt.plot(time, a.f_out)
#plt.title('f_{out}')
#plt.xlabel('Time [s]')
#plt.xlim(xlim1, xlim2)
#plt.ylim(0.99, 1.08)
##plt.tight_layout()
#
#fig12.add_subplot(2,4,6)
#plt.gca().get_yaxis().get_major_formatter().set_useOffset(False) # stop weird exponential notation
#plt.plot(time, a.f_in * a.OEF/p.E_0, label = 'f_{in} * OEF/E_0')
#plt.plot(time, v.HbR/v.CBV * a.f_out, label = 'HbR/CBV * f_{out}')
#plt.title('dHbR/dt components')
#plt.xlabel('Time [s]')
#plt.xlim(xlim1, xlim2)
#plt.ylim(0.99, 1.04)
#plt.legend()
##plt.tight_layout()
#
#fig12.add_subplot(2,4,7)
#plt.gca().get_yaxis().get_major_formatter().set_useOffset(False) # stop weird exponential notation
#plt.plot(time, a.f_in * a.OEF/p.E_0 - v.HbR/v.CBV * a.f_out)
#plt.title('dHbR/dt')
#plt.xlabel('Time [s]')
#plt.xlim(xlim1, xlim2)
#plt.ylim(-0.02, 0.015)
#plt.tight_layout()

#
#fig12.add_subplot(2,4,8)
#plt.gca().get_yaxis().get_major_formatter().set_useOffset(False) # stop weird exponential notation
#plt.plot(time, a.HBR_N)
#plt.title('Normalised HbR')
#plt.xlabel('Time [s]')
#plt.xlim(xlim1, xlim2)
#plt.ylim(0.96, 1.01)
#plt.tight_layout()

#fig2 = plt.figure(dpi=80, edgecolor='k', num = 2) 
#plt.gca().get_yaxis().get_major_formatter().set_useOffset(False) # stop weird exponential notation
#plt.plot(time, a.B_cyt * (a.J_IP3 - a.J_pump + a.J_ER_leak - a.J_TRPV_k/p.r_buff) )
#plt.title('dCak/dt')
#plt.xlim(xlim1, xlim2)
#plt.tight_layout()
#
#fig3 = plt.figure(dpi=80, edgecolor='k', num = 3) 
#plt.gca().get_yaxis().get_major_formatter().set_useOffset(False) # stop weird exponential notation
#plt.plot(time, (p.AA_m * p.AA_max)/(p.AA_m + (v.Ca_k - p.Ca0))**2 * (a.B_cyt * (a.J_IP3 - a.J_pump + a.J_ER_leak - a.J_TRPV_k/p.r_buff) ) + (v.AA_i - v.AA_k)/p.tau_AA)
#plt.title('dAAk/dt')
#plt.xlim(xlim1, xlim2)
#plt.tight_layout()
