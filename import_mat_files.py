# -*- coding: utf-8 -*-
"""
Import the matlab data files
"""

from scipy import io, interpolate
import numpy as np
import parameters as p
import matplotlib.pyplot as plt
from scipy.stats import ttest_ind_from_stats

# Define a class that will contain data (since there are a lot of different arrays in some cases)
class Data():   
    pass

'''Import the Zheng 2010 neural input data'''
# Returns input_data, an array used for InputCase = 'ZhengData', 'ThalamicTrianglesZheng' or 'ZhengFittedParams' for scaling the neuronal stimulus input functions P, Q or T. Based on experimental data from Zheng 2010
# t: time vector
def import_Zheng_neural_data(t):
 
    # If the neuronal stimulation input profile is found from experimental data:
    if p.InputCase == 'ZhengData' or p.InputCase == 'ThalamicTrianglesZheng' or p.InputCase == 'ZhengFittedParams':
        
        # Import the data from mat file
        mat = io.loadmat('./mat_files/Zheng2010_data.mat')
        
        # Extract and flatten data
        neural_tim_vector = mat['neural_tim_vector'].flatten()
        neural_data = mat['neural_data']
   
        sum_neural_wh = [0]*len(neural_tim_vector) # Initialise array
        for animal in range(11):
            for experiment in range(10):
                sum_neural_wh = sum_neural_wh + neural_data[:, p.ISI_idx, p.stim, experiment, animal] # Sum all animal/trial data (110 total)
        mean_neural_wh = sum_neural_wh/110  # Average the neural data
    
        neural_tim_vector_shifted = neural_tim_vector*1e3 + p.startpulse - 20    # Shift so stimulation begins at p.startpulse, -20 so initial spike isn't chopped off
    
        # Interpolate so there is data for all timesteps for NVU
        interp_neural_f = interpolate.interp1d(neural_tim_vector_shifted, mean_neural_wh, fill_value = 0, bounds_error=False)
        interp_neural_wh = interp_neural_f(t)  
    
        if p.double_pulse == 0:  # Remove the second pulse if not wanted
            bool_array = np.in1d(t, p.startpulse + p.lengthpulse+1e3) # find index for where the second pulse begins
            t_idx = np.where(bool_array == True)
            t_idx = t_idx[0][0]
            interp_neural_wh[t_idx:] = 0
    
        # Locus coeruleus pathway - increases input stimulus after some period of time [default = 0]
        if p.LCpathway == 1:
            time_1 = np.linspace(0, p.lengthpulse, 10000) 
            I_LC = 0.0004 * (time_1/1e3)**2 + 1            # Chosen for the shape! Can be modified if needed
            time_1 = time_1 + p.startpulse                 # Shift so starts at p.startpulse
            
            # Interpolate so there is data for all timesteps for NVU
            I_LC_f = interpolate.interp1d(time_1, I_LC, fill_value = 1, bounds_error=False)
            I_LC = I_LC_f(t)  
        else:
             I_LC = 1
        
        input_data = interp_neural_wh * I_LC   # Save input profile for use in P(t), Q(t) or T(t)
                  
    else: # Neuronal stimulation input profile is just a square wave
        input_data = [] # Input data not used
        
    return input_data


'''Import the CBF data from Zheng 2010 experiments'''
# Returns:
# fig - figure of the averaged CBF data
# cbf_time - time variable
# cbf_Zheng - CBF variable
def import_Zheng_cbf_data():
 
    # If the neuronal stimulation input profile is found from experimental data:
    if p.InputCase == 'ZhengData' or p.InputCase == 'ThalamicTrianglesZheng' or p.InputCase == 'ZhengFittedParams':
        
        # Import the data from mat file
        mat = io.loadmat('./mat_files/Zheng2010_data.mat')
        
        # Extract and flatten data
        cbf_tim_vector = mat['cbf_tim_vector'].flatten()
        cbf_data = mat['cbf_data']
       
        sum_cbf = [0]*len(cbf_tim_vector) # Initialise array
        for animal in range(11):
            for experiment in range(10):
                sum_cbf = sum_cbf + cbf_data[:, p.ISI_idx, p.stim, experiment, animal] # Sum all animal/trial data (110 total)
        cbf_Zheng = sum_cbf/110 - 1  # Average the neural data
    
        cbf_time = cbf_tim_vector + p.startpulse*1e-3    # Shift so stimulation begins at p.startpulse, -20 so initial spike isn't chopped off

        fig = plt.figure(figsize=(9,6), dpi=80, edgecolor='k') 
        plt.gca().get_yaxis().get_major_formatter().set_useOffset(False) # stop weird exponential notation
        plt.xlabel('Time [s]')   
        plt.title('CBF (Zheng data)')
        plt.plot(cbf_time,cbf_Zheng)  
       
    else:
        print("Can only import Zheng CBF data if the neuronal stimulation input profile is taken from Zheng data (p.InputCase = 'ZhengData', 'ThalamicTrianglesZheng' or 'ZhengFittedParams')")
        fig, cbf_time, cbf_Zheng = plt.figure(), 0, 0

    return fig, cbf_time, cbf_Zheng

'''Import the data from the Berwick experiments using both air and oxygen'''
# Returns:
# fig1 - figure of the air and oxygen experiments
# d - class containing the time and hemodynamic variables
# area - choose from: 'All', 'Whisker', 'Artery', 'Vein', 'Parenchyma'. Default = 'All'
def import_Berwick_AirOxyData(area='All'):
    
    mat = io.loadmat('./mat_files/Berwick_AirOxy_Data_Full.mat')
    
    # 19 = is the number of sessions (this data is from 5 animals the first four have 4 sessions the last animal has 3 - this is all in order (Animal 1=1:4, Animal 2=5 to 8) - There is one month gap between each of the sessions.
    # 4 = regions from where the responses are from (1 = whisker region, 2 = Artery, 3 = Vein, 4 = Parenchyma)
    # 200 = data length (time)
    # 30 = the number of trials in each of the 19 sessions
       
    # Initialise class to contain the data
    d = Data()
    
    d.time_airOxy = mat['tim'].flatten() - 5                # time vector
    
    '''Find the mean and standard deviation and plot to figure'''
    
    regions = { 'Whisker': 0, 'Artery': 1, 'Vein': 2, 'Parenchyma': 3} # Possible regions
    
    HbOair_array = mat['Hbo_2s_air_Fract']                  # extract the data from the mat file into a 4D array         
    HbRair_array = mat['Hbr_2s_air_Fract']
    HbTair_array = mat['Hbt_2s_air_Fract']
    HbOoxy_array = mat['Hbo_2s_oxy_fract']
    HbRoxy_array = mat['Hbr_2s_oxy_fract']
    HbToxy_array = mat['Hbt_2s_oxy_fract']
    
    if area == 'All':
        HbOair_array = HbOair_array.transpose(3,1,0,2)          # rearrange data so that it can be reshaped in the next step
        HbOair_array = np.reshape(HbOair_array, (19*4*30, 200)) # reshape into a 2D array so that the mean and std can be calculated
        HbRair_array = HbRair_array.transpose(3,1,0,2)
        HbRair_array = np.reshape(HbRair_array, (19*4*30, 200))   
        HbTair_array = HbTair_array.transpose(3,1,0,2)
        HbTair_array = np.reshape(HbTair_array, (19*4*30, 200))   
        HbOoxy_array = HbOoxy_array.transpose(3,1,0,2)
        HbOoxy_array = np.reshape(HbOoxy_array, (19*4*30, 200))         
        HbRoxy_array = HbRoxy_array.transpose(3,1,0,2)
        HbRoxy_array = np.reshape(HbRoxy_array, (19*4*30, 200))   
        HbToxy_array = HbToxy_array.transpose(3,1,0,2)
        HbToxy_array = np.reshape(HbToxy_array, (19*4*30, 200))   
    else:
        HbOair_array = HbOair_array[:,regions[area],:,:].transpose(2,0,1)   # extract data for specific area then rearrange data so that it can be reshaped in the next step
        HbOair_array = np.reshape(HbOair_array, (19*30, 200))               # reshape into a 2D array so that the mean and std can be calculated
        HbRair_array = HbRair_array[:,regions[area],:,:].transpose(2,0,1)   
        HbRair_array = np.reshape(HbRair_array, (19*30, 200))              
        HbTair_array = HbTair_array[:,regions[area],:,:].transpose(2,0,1)   
        HbTair_array = np.reshape(HbTair_array, (19*30, 200))               
        HbOoxy_array = HbOoxy_array[:,regions[area],:,:].transpose(2,0,1)  
        HbOoxy_array = np.reshape(HbOoxy_array, (19*30, 200))               
        HbRoxy_array = HbRoxy_array[:,regions[area],:,:].transpose(2,0,1)   
        HbRoxy_array = np.reshape(HbRoxy_array, (19*30, 200))              
        HbToxy_array = HbToxy_array[:,regions[area],:,:].transpose(2,0,1)  
        HbToxy_array = np.reshape(HbToxy_array, (19*30, 200))               
    
    d.HbOair_mean = np.mean(HbOair_array, axis=0)           # find the mean Hb profile
    d.HbOair_std = np.std(HbOair_array, axis=0)             # find the standard deviation for all the Hb profiles
    d.HbRair_mean = np.mean(HbRair_array, axis=0)
    d.HbRair_std = np.std(HbRair_array, axis=0)
    d.HbTair_mean = np.mean(HbTair_array, axis=0)
    d.HbTair_std = np.std(HbTair_array, axis=0)
    d.HbOoxy_mean = np.mean(HbOoxy_array, axis=0)
    d.HbOoxy_std = np.std(HbOoxy_array, axis=0)
    d.HbRoxy_mean = np.mean(HbRoxy_array, axis=0)
    d.HbRoxy_std = np.std(HbRoxy_array, axis=0)
    d.HbToxy_mean = np.mean(HbToxy_array, axis=0)
    d.HbToxy_std = np.std(HbToxy_array, axis=0)
    
    fig_airOxyComparison = plt.figure(figsize=(12,5), dpi=80, edgecolor='k') 
    fig_airOxyComparison.suptitle("Berwick data ({0:s})".format(area))
    
    # Find ylims based on max/min error bar for either graph
#    ylim1 = min(np.minimum(d.HbRair_mean - d.HbRair_std, d.HbRoxy_mean - d.HbRoxy_std)) - 0.005
#    ylim2 = max(np.maximum(d.HbOair_mean + d.HbOair_std, d.HbOoxy_mean + d.HbOoxy_std)) + 0.005
    ylim1 = 0.91
    ylim2 = 1.12
    
    ax = fig_airOxyComparison.add_subplot(1,2,1)
    ax.get_yaxis().get_major_formatter().set_useOffset(False) # stop weird exponential notation
    ax.set_title('Air')  
    plt.errorbar(d.time_airOxy, d.HbOair_mean, d.HbOair_std, capsize=3, errorevery=12, elinewidth=1, capthick=1, color='r', label = 'HbO')
    plt.errorbar(d.time_airOxy, d.HbRair_mean, d.HbRair_std, capsize=3, errorevery=12, elinewidth=1, capthick=1, color='b', label = 'HbR')
    plt.errorbar(d.time_airOxy, d.HbTair_mean, d.HbTair_std, capsize=3, errorevery=12, elinewidth=1, capthick=1, color='g', label = 'HbT')
    plt.xlabel('Time [s]') 
    plt.ylim(ylim1, ylim2)
    plt.legend()
    
    ax = fig_airOxyComparison.add_subplot(1,2,2)
    ax.get_yaxis().get_major_formatter().set_useOffset(False) # stop weird exponential notation
    ax.set_title('Oxygen')  
    plt.errorbar(d.time_airOxy, d.HbOoxy_mean, d.HbOoxy_std, capsize=3, errorevery=12, elinewidth=1, capthick=1, color='r', label = 'HbO')
    plt.errorbar(d.time_airOxy, d.HbRoxy_mean, d.HbRoxy_std, capsize=3, errorevery=12, elinewidth=1, capthick=1, color='b', label = 'HbR')
    plt.errorbar(d.time_airOxy, d.HbToxy_mean, d.HbToxy_std, capsize=3, errorevery=12, elinewidth=1, capthick=1, color='g', label = 'HbT')
    plt.xlabel('Time [s]') 
    plt.ylim(ylim1, ylim2)
    plt.legend()
    
    '''Find the t and p values at each time step for each Hb and plot to figures'''   
    
    d.HbO_tvalue, d.HbO_pvalue = ttest_ind_from_stats(d.HbOair_mean, d.HbOair_std, 19*4*30, d.HbOoxy_mean, d.HbOoxy_std, 19*4*30)
    d.HbR_tvalue, d.HbR_pvalue = ttest_ind_from_stats(d.HbRair_mean, d.HbRair_std, 19*4*30, d.HbRoxy_mean, d.HbRoxy_std, 19*4*30)
    d.HbT_tvalue, d.HbT_pvalue = ttest_ind_from_stats(d.HbTair_mean, d.HbTair_std, 19*4*30, d.HbToxy_mean, d.HbToxy_std, 19*4*30)    
   
    # Take the absolute value as we don't care about direction
    d.HbO_tvalue = abs(d.HbO_tvalue)
    d.HbR_tvalue = abs(d.HbR_tvalue)
    d.HbT_tvalue = abs(d.HbT_tvalue)
    
    fig_airOxyTTest_HbO = plt.figure(figsize=(10,10), dpi=80, edgecolor='k') 
    fig_airOxyTTest_HbO.suptitle("HbO")
    fig_airOxyTTest_HbO.add_subplot(3,1,1)
    plt.gca().get_yaxis().get_major_formatter().set_useOffset(False) # stop weird exponential notation
    plt.errorbar(d.time_airOxy, d.HbOair_mean, d.HbOair_std, capsize=3, errorevery=12, elinewidth=1.5, capthick=1, label = 'Air')
    plt.errorbar(d.time_airOxy, d.HbOoxy_mean, d.HbOoxy_std, capsize=3, errorevery=12, elinewidth=1.5, capthick=1, label = 'Oxy')
    plt.title('Experimental results')
    plt.xlim(-2,20)  
    plt.legend(loc='best')
    fig_airOxyTTest_HbO.add_subplot(3,1,2)
    plt.gca().get_yaxis().get_major_formatter().set_useOffset(False) # stop weird exponential notation
    plt.plot(d.time_airOxy, d.HbO_tvalue, 'g')
    plt.title('T value')
    plt.xlim(-2,20)  
    fig_airOxyTTest_HbO.add_subplot(3,1,3)
    plt.gca().get_yaxis().get_major_formatter().set_useOffset(False) # stop weird exponential notation
    plt.plot(d.time_airOxy, d.HbO_pvalue, 'purple')
    plt.gca().axhline(y=0.05, label='p = 0.05', c='grey', ls='--') 
    plt.title('P value')
    plt.xlabel('Time [s]') 
    plt.xlim(-2,20)  
    plt.legend()
    
    fig_airOxyTTest_HbR = plt.figure(figsize=(10,10), dpi=80, edgecolor='k') 
    fig_airOxyTTest_HbR.suptitle("HbR")
    fig_airOxyTTest_HbR.add_subplot(3,1,1)
    plt.gca().get_yaxis().get_major_formatter().set_useOffset(False) # stop weird exponential notation
    plt.errorbar(d.time_airOxy, d.HbRair_mean, d.HbRair_std, capsize=3, errorevery=12, elinewidth=1.5, capthick=1, label = 'Air')
    plt.errorbar(d.time_airOxy, d.HbRoxy_mean, d.HbRoxy_std, capsize=3, errorevery=12, elinewidth=1.5, capthick=1, label = 'Oxy')
    plt.title('Experimental results')
    plt.xlim(-2,20)  
    plt.legend(loc='best')
    fig_airOxyTTest_HbR.add_subplot(3,1,2)
    plt.gca().get_yaxis().get_major_formatter().set_useOffset(False) # stop weird exponential notation
    plt.plot(d.time_airOxy, d.HbR_tvalue, 'g')
    plt.title('T value')
    plt.xlim(-2,20)  
    fig_airOxyTTest_HbR.add_subplot(3,1,3)
    plt.gca().get_yaxis().get_major_formatter().set_useOffset(False) # stop weird exponential notation
    plt.plot(d.time_airOxy, d.HbR_pvalue, 'purple')
    plt.gca().axhline(y=0.05, label='p = 0.05', c='grey', ls='--') 
    plt.title('P value')
    plt.xlabel('Time [s]') 
    plt.xlim(-2,20)  
    plt.legend()
    
    fig_airOxyTTest_HbT = plt.figure(figsize=(10,10), dpi=80, edgecolor='k') 
    fig_airOxyTTest_HbT.suptitle("HbT")
    fig_airOxyTTest_HbT.add_subplot(3,1,1)
    plt.gca().get_yaxis().get_major_formatter().set_useOffset(False) # stop weird exponential notation
    plt.errorbar(d.time_airOxy, d.HbTair_mean, d.HbTair_std, capsize=3, errorevery=12, elinewidth=1.5, capthick=1, label = 'Air')
    plt.errorbar(d.time_airOxy, d.HbToxy_mean, d.HbToxy_std, capsize=3, errorevery=12, elinewidth=1.5, capthick=1, label = 'Oxy')
    plt.title('Experimental results')
    plt.xlim(-2,20)  
    plt.legend(loc='best')
    fig_airOxyTTest_HbT.add_subplot(3,1,2)
    plt.gca().get_yaxis().get_major_formatter().set_useOffset(False) # stop weird exponential notation
    plt.plot(d.time_airOxy, d.HbT_tvalue, 'g')
    plt.title('T value')
    plt.xlim(-2,20)  
    fig_airOxyTTest_HbT.add_subplot(3,1,3)
    plt.gca().get_yaxis().get_major_formatter().set_useOffset(False) # stop weird exponential notation
    plt.plot(d.time_airOxy, d.HbT_pvalue, 'purple')
    plt.gca().axhline(y=0.05, label='p = 0.05', c='grey', ls='--')
    plt.title('P value')
    plt.xlabel('Time [s]') 
    plt.xlim(-2,20) 
    plt.legend()
    
    return fig_airOxyComparison, fig_airOxyTTest_HbO, fig_airOxyTTest_HbR, fig_airOxyTTest_HbT, d
      
'''Import the LNAME/7NI data from Berwick et al'''
# Returns:
# fig1 - figure of the whisker experiments
# fig2 - figure of the optogenetic experiments
# d - class containing the time and hemodynamic variables
# area - choose from: 'All', 'Whisker', 'Artery', 'Vein', 'Parenchyma'. Default = 'All'
def import_Berwick_LNAME_7NI_data(area='All'):
    
    mat = io.loadmat('./mat_files/LNAME_data_preliminary.mat')
    
    # 8 = is the number of animals.
    # 4 = regions from where the responses are from (1 = whisker region, 2 = Artery, 3 = Vein, 4 = Parenchyma)
    # 200 = data length (time)
    
    # Initialise class to contain the data
    d = Data()
    
    d.time_LNAME = mat['tim'].flatten() - 5  # Time array
    
    regions = { 'Whisker': 0, 'Artery': 1, 'Vein': 2, 'Parenchyma': 3} # Possible regions
    saturation = { 'Whisker': 70, 'Artery': 80, 'Vein': 60, 'Parenchyma': 70} # saturation for each region, used for scaling the micromolar changes into fractional changes

    # Initialise arrays to contain the data for the 4 regions
    HbO_wh_preLNAME_all = [0]*4
    HbO_wh_postLNAME_all = [0]*4
    HbO_op_preLNAME_all = [0]*4
    HbO_op_postLNAME_all = [0]*4
    HbR_wh_preLNAME_all = [0]*4
    HbR_wh_postLNAME_all = [0]*4
    HbR_op_preLNAME_all = [0]*4
    HbR_op_postLNAME_all = [0]*4 
    HbT_wh_preLNAME_all = [0]*4
    HbT_wh_postLNAME_all = [0]*4
    HbT_op_preLNAME_all = [0]*4
    HbT_op_postLNAME_all = [0]*4
    
    # Extract and scale data for each region and each experiment type
    for r in regions:
        HbO_wh_preLNAME_all[regions[r]] = mat['hbo_whisker_pre'][:,regions[r],:]/saturation[r] + 1 
        HbO_wh_postLNAME_all[regions[r]] = mat['hbo_whisker_post'][:,regions[r],:]/saturation[r] + 1   
        HbO_op_preLNAME_all[regions[r]] = mat['hbo_opto_pre'][:,regions[r],:]/saturation[r] + 1   
        HbO_op_postLNAME_all[regions[r]] = mat['hbo_opto_post'][:,regions[r],:]/saturation[r] + 1     
        HbR_wh_preLNAME_all[regions[r]] = mat['hbr_whisker_pre'][:,regions[r],:]/(100-saturation[r]) + 1   
        HbR_wh_postLNAME_all[regions[r]] = mat['hbr_whisker_post'][:,regions[r],:]/(100-saturation[r]) + 1       
        HbR_op_preLNAME_all[regions[r]] = mat['hbr_opto_pre'][:,regions[r],:]/(100-saturation[r]) + 1   
        HbR_op_postLNAME_all[regions[r]] = mat['hbr_opto_post'][:,regions[r],:]/(100-saturation[r]) + 1   
        HbT_wh_preLNAME_all[regions[r]] = mat['hbt_whisker_pre'][:,regions[r],:]/100 + 1   
        HbT_wh_postLNAME_all[regions[r]] = mat['hbt_whisker_post'][:,regions[r],:]/100 + 1       
        HbT_op_preLNAME_all[regions[r]] = mat['hbt_opto_pre'][:,regions[r],:]/100 + 1   
        HbT_op_postLNAME_all[regions[r]] = mat['hbt_opto_post'][:,regions[r],:]/100 + 1   
      
    
    # Extract and shape the arrays for each case
    if area == 'All':   # stack data from all regions into a single 2D array 
        HbO_wh_preLNAME = np.vstack((HbO_wh_preLNAME_all[0], HbO_wh_preLNAME_all[1], HbO_wh_preLNAME_all[2], HbO_wh_preLNAME_all[3]))
        HbO_wh_postLNAME =np.vstack((HbO_wh_postLNAME_all[0], HbO_wh_postLNAME_all[1], HbO_wh_postLNAME_all[2], HbO_wh_postLNAME_all[3]))
        HbO_op_preLNAME = np.vstack((HbO_op_preLNAME_all[0], HbO_op_preLNAME_all[1], HbO_op_preLNAME_all[2], HbO_op_preLNAME_all[3]))
        HbO_op_postLNAME = np.vstack((HbO_op_postLNAME_all[0], HbO_op_postLNAME_all[1], HbO_op_postLNAME_all[2], HbO_op_postLNAME_all[3]))
        HbR_wh_preLNAME = np.vstack((HbR_wh_preLNAME_all[0], HbR_wh_preLNAME_all[1], HbR_wh_preLNAME_all[2], HbR_wh_preLNAME_all[3]))
        HbR_wh_postLNAME = np.vstack((HbR_wh_postLNAME_all[0], HbR_wh_postLNAME_all[1], HbR_wh_postLNAME_all[2], HbR_wh_postLNAME_all[3]))
        HbR_op_preLNAME = np.vstack((HbR_op_preLNAME_all[0], HbR_op_preLNAME_all[1], HbR_op_preLNAME_all[2], HbR_op_preLNAME_all[3]))
        HbR_op_postLNAME = np.vstack((HbR_op_postLNAME_all[0], HbR_op_postLNAME_all[1], HbR_op_postLNAME_all[2], HbR_op_postLNAME_all[3]))
        HbT_wh_preLNAME = np.vstack((HbT_wh_preLNAME_all[0], HbT_wh_preLNAME_all[1], HbT_wh_preLNAME_all[2], HbT_wh_preLNAME_all[3]))
        HbT_wh_postLNAME = np.vstack((HbT_wh_postLNAME_all[0], HbT_wh_postLNAME_all[1], HbT_wh_postLNAME_all[2], HbT_wh_postLNAME_all[3]))
        HbT_op_preLNAME = np.vstack((HbT_op_preLNAME_all[0], HbT_op_preLNAME_all[1], HbT_op_preLNAME_all[2], HbT_op_preLNAME_all[3]))
        HbT_op_postLNAME = np.vstack((HbT_op_postLNAME_all[0], HbT_op_postLNAME_all[1], HbT_op_postLNAME_all[2], HbT_op_postLNAME_all[3]))
    else:               # extract data for specific area 
        HbO_wh_preLNAME = HbO_wh_preLNAME_all[regions[area]]
        HbO_wh_postLNAME = HbO_wh_postLNAME_all[regions[area]]  
        HbO_op_preLNAME = HbO_op_preLNAME_all[regions[area]]
        HbO_op_postLNAME = HbO_op_postLNAME_all[regions[area]] 
        HbR_wh_preLNAME = HbR_wh_preLNAME_all[regions[area]]
        HbR_wh_postLNAME = HbR_wh_postLNAME_all[regions[area]]     
        HbR_op_preLNAME = HbR_op_preLNAME_all[regions[area]] 
        HbR_op_postLNAME = HbR_op_postLNAME_all[regions[area]]  
        HbT_wh_preLNAME = HbT_wh_preLNAME_all[regions[area]] 
        HbT_wh_postLNAME = HbT_wh_postLNAME_all[regions[area]]  
        HbT_op_preLNAME = HbT_op_preLNAME_all[regions[area]] 
        HbT_op_postLNAME = HbT_op_postLNAME_all[regions[area]]
    
    # Find the means and standard deviations
    d.HbO_wh_preLNAME_mean = np.mean(HbO_wh_preLNAME, axis=0)          
    d.HbO_wh_preLNAME_std = np.std(HbO_wh_preLNAME, axis=0)             
    d.HbO_wh_postLNAME_mean = np.mean(HbO_wh_postLNAME, axis=0)          
    d.HbO_wh_postLNAME_std = np.std(HbO_wh_postLNAME, axis=0)   
    d.HbO_op_preLNAME_mean = np.mean(HbO_op_preLNAME, axis=0)          
    d.HbO_op_preLNAME_std = np.std(HbO_op_preLNAME, axis=0)   
    d.HbO_op_postLNAME_mean = np.mean(HbO_op_postLNAME, axis=0)          
    d.HbO_op_postLNAME_std = np.std(HbO_op_postLNAME, axis=0) 
    
    d.HbR_wh_preLNAME_mean = np.mean(HbR_wh_preLNAME, axis=0)          
    d.HbR_wh_preLNAME_std = np.std(HbR_wh_preLNAME, axis=0)             
    d.HbR_wh_postLNAME_mean = np.mean(HbR_wh_postLNAME, axis=0)          
    d.HbR_wh_postLNAME_std = np.std(HbR_wh_postLNAME, axis=0)   
    d.HbR_op_preLNAME_mean = np.mean(HbR_op_preLNAME, axis=0)          
    d.HbR_op_preLNAME_std = np.std(HbR_op_preLNAME, axis=0)   
    d.HbR_op_postLNAME_mean = np.mean(HbR_op_postLNAME, axis=0)          
    d.HbR_op_postLNAME_std = np.std(HbR_op_postLNAME, axis=0)  
    
    d.HbT_wh_preLNAME_mean = np.mean(HbT_wh_preLNAME, axis=0)          
    d.HbT_wh_preLNAME_std = np.std(HbT_wh_preLNAME, axis=0)             
    d.HbT_wh_postLNAME_mean = np.mean(HbT_wh_postLNAME, axis=0)          
    d.HbT_wh_postLNAME_std = np.std(HbT_wh_postLNAME, axis=0)   
    d.HbT_op_preLNAME_mean = np.mean(HbT_op_preLNAME, axis=0)          
    d.HbT_op_preLNAME_std = np.std(HbT_op_preLNAME, axis=0)   
    d.HbT_op_postLNAME_mean = np.mean(HbT_op_postLNAME, axis=0)          
    d.HbT_op_postLNAME_std = np.std(HbT_op_postLNAME, axis=0)  

    fig = plt.figure(figsize=(12,10), dpi=80, edgecolor='k') 
    fig.suptitle("LNAME data ({0:s})".format(area))
    
    # Find ylims based on max/min error bar for either graph
    ylim1 = min(min(d.HbR_wh_preLNAME_mean - d.HbR_wh_preLNAME_std), min(d.HbR_wh_postLNAME_mean - d.HbR_wh_postLNAME_std), min(d.HbR_op_preLNAME_mean - d.HbR_op_preLNAME_std), min(d.HbR_op_postLNAME_mean - d.HbR_op_postLNAME_std)) - 0.002
    ylim2 = max(max(d.HbO_wh_preLNAME_mean + d.HbO_wh_preLNAME_std), max(d.HbO_wh_postLNAME_mean + d.HbO_wh_postLNAME_std), max(d.HbO_op_preLNAME_mean + d.HbO_op_preLNAME_std), max(d.HbO_op_postLNAME_mean + d.HbO_op_postLNAME_std)) + 0.002
    
    ax = fig.add_subplot(2,2,1)
    ax.get_yaxis().get_major_formatter().set_useOffset(False) # stop weird exponential notation
    ax.set_title('Whisker Response - Pre LNAME')  
    plt.errorbar(d.time_LNAME, d.HbO_wh_preLNAME_mean, d.HbO_wh_preLNAME_std, capsize=3, errorevery=12, elinewidth=1, capthick=1, color='r', label = 'HbO')
    plt.errorbar(d.time_LNAME, d.HbR_wh_preLNAME_mean, d.HbR_wh_preLNAME_std, capsize=3, errorevery=12, elinewidth=1, capthick=1, color='b', label = 'HbR')
    plt.errorbar(d.time_LNAME, d.HbT_wh_preLNAME_mean, d.HbT_wh_preLNAME_std, capsize=3, errorevery=12, elinewidth=1, capthick=1, color='g', label = 'HbT')
    plt.xlabel('Time [s]') 
    plt.ylim(ylim1, ylim2)
    plt.legend()
    
    ax = fig.add_subplot(2,2,2)
    ax.get_yaxis().get_major_formatter().set_useOffset(False) # stop weird exponential notation
    ax.set_title('Whisker Response - Post LNAME')  
    plt.errorbar(d.time_LNAME, d.HbO_wh_postLNAME_mean, d.HbO_wh_postLNAME_std, capsize=3, errorevery=12, elinewidth=1, capthick=1, color='r', label = 'HbO')
    plt.errorbar(d.time_LNAME, d.HbR_wh_postLNAME_mean, d.HbR_wh_postLNAME_std, capsize=3, errorevery=12, elinewidth=1, capthick=1, color='b', label = 'HbR')
    plt.errorbar(d.time_LNAME, d.HbT_wh_postLNAME_mean, d.HbT_wh_postLNAME_std, capsize=3, errorevery=12, elinewidth=1, capthick=1, color='g', label = 'HbT')
    plt.xlabel('Time [s]') 
    plt.ylim(ylim1, ylim2)
    plt.legend()   
    
    ax = fig.add_subplot(2,2,3)
    ax.get_yaxis().get_major_formatter().set_useOffset(False) # stop weird exponential notation
    ax.set_title('Opto Response - Pre LNAME')  
    plt.errorbar(d.time_LNAME, d.HbO_op_preLNAME_mean, d.HbO_op_preLNAME_std, capsize=3, errorevery=12, elinewidth=1, capthick=1, color='r', label = 'HbO')
    plt.errorbar(d.time_LNAME, d.HbR_op_preLNAME_mean, d.HbR_op_preLNAME_std, capsize=3, errorevery=12, elinewidth=1, capthick=1, color='b', label = 'HbR')
    plt.errorbar(d.time_LNAME, d.HbT_op_preLNAME_mean, d.HbT_op_preLNAME_std, capsize=3, errorevery=12, elinewidth=1, capthick=1, color='g', label = 'HbT')
    plt.xlabel('Time [s]') 
    plt.ylim(ylim1, ylim2)
    plt.legend()  

    ax = fig.add_subplot(2,2,4)
    ax.get_yaxis().get_major_formatter().set_useOffset(False) # stop weird exponential notation
    ax.set_title('Opto Response - Post LNAME')  
    plt.errorbar(d.time_LNAME, d.HbO_op_postLNAME_mean, d.HbO_op_postLNAME_std, capsize=3, errorevery=12, elinewidth=1, capthick=1, color='r', label = 'HbO')
    plt.errorbar(d.time_LNAME, d.HbR_op_postLNAME_mean, d.HbR_op_postLNAME_std, capsize=3, errorevery=12, elinewidth=1, capthick=1, color='b', label = 'HbR')
    plt.errorbar(d.time_LNAME, d.HbT_op_postLNAME_mean, d.HbT_op_postLNAME_std, capsize=3, errorevery=12, elinewidth=1, capthick=1, color='g', label = 'HbT')
    plt.xlabel('Time [s]') 
    plt.ylim(ylim1, ylim2)
    plt.legend()  

    ## Load, extract, flatten and scale the data
#    mat = io.loadmat('./mat_files/Hb_whisker_pre7NI.mat')
#    d.HbO_wh_pre7ni = mat['HbOPreMicromolarchange'].flatten()/70 + 1
#    d.HbR_wh_pre7ni = mat['HbRPreMicromolarchange'].flatten()/30 + 1
#    d.HbT_wh_pre7ni = mat['HbTPremicromolarchange'].flatten()/100 + 1  
#    d.time_HbO_wh_pre7ni = mat['Time6'].flatten() - 5
#    d.time_HbR_wh_pre7ni = mat['Time7'].flatten() - 5
#    d.time_HbT_wh_pre7ni = mat['Time8'].flatten() - 5    
#
#    mat = io.loadmat('./mat_files/Hb_whisker_post7NI.mat')  
#    d.HbO_wh_post7ni = mat['HbOPost7NIMicromolarchange'].flatten()/70 + 1
#    d.HbR_wh_post7ni = mat['HbRPost7NIMicromolarchange'].flatten()/30 + 1
#    d.HbT_wh_post7ni = mat['HbTPost7NIMicromolarchange'].flatten()/100 + 1  
#    d.time_HbO_wh_post7ni = mat['TIme9'].flatten() - 5
#    d.time_HbR_wh_post7ni = mat['Time10'].flatten() - 5
#    d.time_HbT_wh_post7ni = mat['Time11'].flatten() - 5
#    
#    mat = io.loadmat('./mat_files/Hb_opto_pre7NI.mat')  
#    d.HbO_op_pre7ni = mat['HbO'].flatten()/70 + 1
#    d.HbR_op_pre7ni = mat['HbR'].flatten()/30 + 1
#    d.HbT_op_pre7ni = mat['HbT'].flatten()/100 + 1  
#    d.time_HbO_op_pre7ni = mat['TimeHbO'].flatten() - 5
#    d.time_HbR_op_pre7ni = mat['TimeHbR'].flatten() - 5
#    d.time_HbT_op_pre7ni = mat['TimeHbT'].flatten() - 5
#
#    mat = io.loadmat('./mat_files/Hb_opto_post7NI.mat')   
#    d.HbO_op_post7ni = mat['HbOMicromolarChange'].flatten()/70 + 1
#    d.HbR_op_post7ni = mat['HbRMicromolarChange'].flatten()/30 + 1
#    d.HbT_op_post7ni = mat['HbTMicromolarChange'].flatten()/100 + 1  
#    d.time_HbO_op_post7ni = mat['TimeHbO'].flatten() - 5
#    d.time_HbR_op_post7ni = mat['TimeHbR'].flatten() - 5
#    d.time_HbT_op_post7ni = mat['TimeHbT'].flatten() - 5
#
#    mat = io.loadmat('./mat_files/Hb_whisker_preLNAME.mat')   
#    d.HbO_wh_preLNAME = mat['HbO'].flatten()/70 + 1
#    d.HbR_wh_preLNAME= mat['HbR'].flatten()/30 + 1
#    d.HbT_wh_preLNAME = mat['HbT'].flatten()/100 + 1  
#    d.time_HbO_wh_preLNAME = mat['TimeHbO'].flatten() - 4.5
#    d.time_HbR_wh_preLNAME = mat['TimeHbR'].flatten() - 4.5
#    d.time_HbT_wh_preLNAME = mat['TimeHbT'].flatten() - 4.5
#
#    mat = io.loadmat('./mat_files/Hb_whisker_postLNAME.mat')   
#    d.HbO_wh_postLNAME = mat['HbOPostMicromolarchange'].flatten()/70 + 1
#    d.HbR_wh_postLNAME= mat['HbRPostMicromolarchange'].flatten()/30 + 1
#    d.HbT_wh_postLNAME = mat['HbTPostMicromolarchange'].flatten()/100 + 1  
#    d.time_HbO_wh_postLNAME = mat['Time3'].flatten() - 4.5
#    d.time_HbR_wh_postLNAME = mat['Time4'].flatten() - 4.5
#    d.time_HbT_wh_postLNAME = mat['Time5'].flatten() - 4.5
#
#    mat = io.loadmat('./mat_files/Hb_opto_preLNAME.mat')   
#    d.HbO_op_preLNAME = mat['HbOMicromolarchange'].flatten()/70 + 1
#    d.HbR_op_preLNAME= mat['HbRMicromolarchange'].flatten()/30 + 1
#    d.HbT_op_preLNAME = mat['HbTMicromolarchange'].flatten()/100 + 1  
#    d.time_HbO_op_preLNAME = mat['Timeopt1'].flatten() - 4.5
#    d.time_HbR_op_preLNAME = mat['Timeopt2'].flatten() - 4.5
#    d.time_HbT_op_preLNAME = mat['TImeopt3'].flatten() - 4.5
#
#    mat = io.loadmat('./mat_files/Hb_opto_postLNAME.mat')
#    d.HbO_op_postLNAME = mat['HbOMicromolarchange'].flatten()/70 + 1
#    d.HbR_op_postLNAME= mat['HbRMicromolarchange'].flatten()/30 + 1
#    d.HbT_op_postLNAME = mat['HbTMicromolarchange'].flatten()/100 + 1  
#    d.time_HbO_op_postLNAME = mat['Timeopt4'].flatten() - 4.5
#    d.time_HbR_op_postLNAME = mat['Timeopt5'].flatten() - 4.5
#    d.time_HbT_op_postLNAME = mat['Timeopt6'].flatten() - 4.5
#    
#    ## Plot the whisker results
#    fig1 = plt.figure(figsize=(10,6), dpi=80, edgecolor='k') 
#    fig1.suptitle("Berwick data: Whisker response")
#
#    ax = fig1.add_subplot(2,2,1)
#    ax.get_yaxis().get_major_formatter().set_useOffset(False) # stop weird exponential notation
#    ax.set_title('Pre 7NI')
#    ax.plot(d.time_HbO_wh_pre7ni,d.HbO_wh_pre7ni, 'r', label = 'HbO')
#    ax.plot(d.time_HbR_wh_pre7ni,d.HbR_wh_pre7ni, 'b', label = 'HbR')
#    ax.plot(d.time_HbT_wh_pre7ni,d.HbT_wh_pre7ni, 'g', label = 'HbT')
#    plt.xlabel('Time [s]') 
#    plt.legend()
#    
#    ax = fig1.add_subplot(2,2,2)
#    ax.get_yaxis().get_major_formatter().set_useOffset(False) # stop weird exponential notation
#    ax.set_title('Post 7NI')
#    ax.plot(d.time_HbO_wh_post7ni,d.HbO_wh_post7ni, 'r', label = 'HbO')
#    ax.plot(d.time_HbR_wh_post7ni,d.HbR_wh_post7ni, 'b', label = 'HbR')
#    ax.plot(d.time_HbT_wh_post7ni,d.HbT_wh_post7ni, 'g', label = 'HbT')
#    plt.xlabel('Time [s]') 
#    plt.legend()
#
#    ax = fig1.add_subplot(2,2,3)
#    ax.get_yaxis().get_major_formatter().set_useOffset(False) # stop weird exponential notation
#    ax.set_title('Pre L-NAME')
#    ax.plot(d.time_HbO_wh_preLNAME,d.HbO_wh_preLNAME, 'r', label = 'HbO')
#    ax.plot(d.time_HbR_wh_preLNAME,d.HbR_wh_preLNAME, 'b', label = 'HbR')
#    ax.plot(d.time_HbT_wh_preLNAME,d.HbT_wh_preLNAME, 'g', label = 'HbT')
#    plt.xlabel('Time [s]') 
#    plt.legend()
#
#    ax = fig1.add_subplot(2,2,4)
#    ax.get_yaxis().get_major_formatter().set_useOffset(False) # stop weird exponential notation
#    ax.set_title('Post L-NAME')
#    ax.plot(d.time_HbO_wh_postLNAME,d.HbO_wh_postLNAME, 'r', label = 'HbO')
#    ax.plot(d.time_HbR_wh_postLNAME,d.HbR_wh_postLNAME, 'b', label = 'HbR')
#    ax.plot(d.time_HbT_wh_postLNAME,d.HbT_wh_postLNAME, 'g', label = 'HbT')
#    plt.xlabel('Time [s]') 
#    plt.legend()    
#    
#    plt.tight_layout()
#    
#    # Plot the optogenetic results
#    fig2 = plt.figure(figsize=(10,6), dpi=80, edgecolor='k') 
#    fig2.suptitle("Berwick data: Opto response")
#
#    ax = fig2.add_subplot(2,2,1)
#    ax.get_yaxis().get_major_formatter().set_useOffset(False) # stop weird exponential notation
#    ax.set_title('Pre 7NI')
#    ax.plot(d.time_HbO_op_pre7ni,d.HbO_op_pre7ni, 'r', label = 'HbO')
#    ax.plot(d.time_HbR_op_pre7ni,d.HbR_op_pre7ni, 'b', label = 'HbR')
#    ax.plot(d.time_HbT_op_pre7ni,d.HbT_op_pre7ni, 'g', label = 'HbT')
#    plt.xlabel('Time [s]') 
#    plt.legend()
#    
#    ax = fig2.add_subplot(2,2,2)
#    ax.get_yaxis().get_major_formatter().set_useOffset(False) # stop weird exponential notation
#    ax.set_title('Post 7NI')
#    ax.plot(d.time_HbO_op_post7ni,d.HbO_op_post7ni, 'r', label = 'HbO')
#    ax.plot(d.time_HbR_op_post7ni,d.HbR_op_post7ni, 'b', label = 'HbR')
#    ax.plot(d.time_HbT_op_post7ni,d.HbT_op_post7ni, 'g', label = 'HbT')
#    plt.xlabel('Time [s]') 
#    plt.legend()
#
#    ax = fig2.add_subplot(2,2,3)
#    ax.get_yaxis().get_major_formatter().set_useOffset(False) # stop weird exponential notation
#    ax.set_title('Pre L-NAME')
#    ax.plot(d.time_HbO_op_preLNAME,d.HbO_op_preLNAME, 'r', label = 'HbO')
#    ax.plot(d.time_HbR_op_preLNAME,d.HbR_op_preLNAME, 'b', label = 'HbR')
#    ax.plot(d.time_HbT_op_preLNAME,d.HbT_op_preLNAME, 'g', label = 'HbT')
#    plt.xlabel('Time [s]') 
#    plt.legend()
#
#    ax = fig2.add_subplot(2,2,4)
#    ax.get_yaxis().get_major_formatter().set_useOffset(False) # stop weird exponential notation
#    ax.set_title('Post L-NAME')
#    ax.plot(d.time_HbO_op_postLNAME,d.HbO_op_postLNAME, 'r', label = 'HbO')
#    ax.plot(d.time_HbR_op_postLNAME,d.HbR_op_postLNAME, 'b', label = 'HbR')
#    ax.plot(d.time_HbT_op_postLNAME,d.HbT_op_postLNAME, 'g', label = 'HbT')
#    plt.xlabel('Time [s]') 
#    plt.legend()    
#    
#    plt.tight_layout()
    
    return fig, d

'''Extract the experimental data from the LNAME/7NI data that corresponds to the current model simulation'''
# TODO: change to match new data
# Returns the 3 time variables and 3 hemodynamic variables
# d: class containing the experimental data 
def set_corresponding_LNAME_7NI_experimental_data(d):
    
    time_LNAME = d.time_LNAME 
    
    if p.NeuronType == 'whisker' and p.NOswitch == 'normal':    
        exp_HbO = d.HbO_wh_preLNAME_mean
        exp_HbR = d.HbR_wh_preLNAME_mean
        exp_HbT = d.HbT_wh_preLNAME_mean  
        exp_HbO_std = d.HbO_wh_preLNAME_std
        exp_HbR_std = d.HbR_wh_preLNAME_std 
        exp_HbT_std = d.HbT_wh_preLNAME_std 
    elif p.NeuronType == 'whisker' and p.NOswitch == 'LNAME': 
        exp_HbO = d.HbO_wh_postLNAME_mean
        exp_HbR = d.HbR_wh_postLNAME_mean
        exp_HbT = d.HbT_wh_postLNAME_mean 
        exp_HbO_std = d.HbO_wh_postLNAME_std
        exp_HbR_std = d.HbR_wh_postLNAME_std
        exp_HbT_std = d.HbT_wh_postLNAME_std  
    elif p.NeuronType == 'opto' and p.NOswitch == 'normal':    
        exp_HbO = d.HbO_op_preLNAME_mean
        exp_HbR = d.HbR_op_preLNAME_mean
        exp_HbT = d.HbT_op_preLNAME_mean
        exp_HbO_std = d.HbO_op_preLNAME_std
        exp_HbR_std = d.HbR_op_preLNAME_std
        exp_HbT_std = d.HbT_op_preLNAME_std 
    elif p.NeuronType == 'opto' and p.NOswitch == 'LNAME':
        exp_HbO = d.HbO_op_postLNAME_mean
        exp_HbR = d.HbR_op_postLNAME_mean
        exp_HbT = d.HbT_op_postLNAME_mean  
        exp_HbO_std = d.HbO_op_postLNAME_std
        exp_HbR_std = d.HbR_op_postLNAME_std
        exp_HbT_std = d.HbT_op_postLNAME_std      
#    if p.NeuronType == 'whisker' and p.NOswitch == 'normal':
#        time_HbO = d.time_HbO_wh_preLNAME
#        time_HbR = d.time_HbR_wh_preLNAME
#        time_HbT = d.time_HbT_wh_preLNAME   
#        exp_HbO = d.HbO_wh_preLNAME
#        exp_HbR = d.HbR_wh_preLNAME
#        exp_HbT = d.HbT_wh_preLNAME  
#    elif p.NeuronType == 'whisker' and p.NOswitch == '7NI':
#        time_HbO = d.time_HbO_wh_post7ni
#        time_HbR = d.time_HbR_wh_post7ni
#        time_HbT = d.time_HbT_wh_post7ni   
#        exp_HbO = d.HbO_wh_post7ni
#        exp_HbR = d.HbR_wh_post7ni
#        exp_HbT = d.HbT_wh_post7ni
#    elif p.NeuronType == 'whisker' and p.NOswitch == 'LNAME':
#        time_HbO = d.time_HbO_wh_postLNAME
#        time_HbR = d.time_HbR_wh_postLNAME
#        time_HbT = d.time_HbT_wh_postLNAME   
#        exp_HbO = d.HbO_wh_postLNAME
#        exp_HbR = d.HbR_wh_postLNAME
#        exp_HbT = d.HbT_wh_postLNAME    
#    elif p.NeuronType == 'opto' and p.NOswitch == 'normal':
#        time_HbO = d.time_HbO_op_preLNAME
#        time_HbR = d.time_HbR_op_preLNAME
#        time_HbT = d.time_HbT_op_preLNAME   
#        exp_HbO = d.HbO_op_preLNAME
#        exp_HbR = d.HbR_op_preLNAME
#        exp_HbT = d.HbT_op_preLNAME  
#    elif p.NeuronType == 'opto' and p.NOswitch == '7NI':
#        time_HbO = d.time_HbO_op_post7ni
#        time_HbR = d.time_HbR_op_post7ni
#        time_HbT = d.time_HbT_op_post7ni   
#        exp_HbO = d.HbO_op_post7ni
#        exp_HbR = d.HbR_op_post7ni
#        exp_HbT = d.HbT_op_post7ni
#    elif p.NeuronType == 'opto' and p.NOswitch == 'LNAME':
#        time_HbO = d.time_HbO_op_postLNAME
#        time_HbR = d.time_HbR_op_postLNAME
#        time_HbT = d.time_HbT_op_postLNAME   
#        exp_HbO = d.HbO_op_postLNAME
#        exp_HbR = d.HbR_op_postLNAME
#        exp_HbT = d.HbT_op_postLNAME  
        
    return time_LNAME, exp_HbO, exp_HbR, exp_HbT, exp_HbO_std, exp_HbR_std, exp_HbT_std








