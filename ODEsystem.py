# -*- coding: utf-8 -*-
"""
The differential equations of the NVU model
"""

from algebraic_variables import set_algebraic_variables
from state_variable_names import initialise_var_names
import parameters as p
import time

# Returns du, the RHS array of the differential equation set for use in the solver ODEINT
# u: array of state variables 
# t: time vector
# idx: dictionary containing the state variable names ('keys') and their corresponding index ('values')
# input_data: experimental data used for the neuronal input functions
# total_waveform: an array containing the fast triangular input pulses when using the thalamic input T(t)
# start_time: start time of simulation, used to calculate the time since beginning
def func(u ,t, idx, input_data, total_waveform, start_time):
    
    # Initialise variable names for clarity 
    v = initialise_var_names(u, idx, 'in function')
    
    # Set the algebraic_variables
    a = set_algebraic_variables(u, t, idx, input_data, total_waveform, 'in function')
    
    # RHS of ODEs
    du = [None]*len(idx) # Initialise with correct size
       
    if p.InputCase == 'ThalamicTriangles' or p.InputCase == 'ThalamicTrianglesZheng' or p.InputCase == 'ThalamicSquarePulse': # Using thalamic input T(t) rather than P(t), Q(t)
        du[idx['E_t']] =  1 / p.tau_e * (-v.E_t + a.f_arg_e)            
        du[idx['I_t']] =  1 / p.tau_i * (-v.I_t + a.f_arg_i)          
    else:
        du[idx['E_t']] =  1 / p.tau_e * (-v.E_t + (a.k_e - p.r_e * v.E_t) * a.s_arg_e)                 
        du[idx['I_t']] =  1 / p.tau_i * (-v.I_t + (a.k_i - p.r_i * v.I_t) * a.s_arg_i)             
    
    du[idx['K_e']] = - p.beta_K_e * (v.K_e - p.K_eBase) + p.alpha_K_e * p.beta_K_e * ((abs(v.E_t - v.I_t) - p.EImin) / (p.EI_relative - p.EImin))  
    du[idx['Na_sa']] = - p.beta_Na_sa * (v.Na_sa - p.Na_saBase) + p.alpha_Na_sa * p.beta_Na_sa * ((abs(v.E_t - v.I_t) - p.EImin)/ (p.EI_relative - p.EImin))                   
    du[idx['Na_d']] =  - p.beta_Na_d * (v.Na_d - p.Na_dBase) +  p.alpha_Na_d * p.beta_Na_d * ((abs(v.E_t - v.I_t) - p.EImin)/ (p.EI_relative - p.EImin))
    du[idx['O2']] = a.J_O2_vascular - a.J_O2_background - a.J_O2_pump            
    du[idx['CBV']] = 1/(p.tau_MTT + p.tau_TAT) * ( a.f_in  - v.CBV**(1/p.d) )          
    du[idx['HbR']] = 1/p.tau_MTT * ( a.f_in * a.OEF / p.E_0 - v.HbR/v.CBV * a.f_out )  
    
    du[idx['Ca_n']] = (a.I_Ca_tot / (2 * p.Farad * p.V_spine) - (p.k_ex * (v.Ca_n - p.Ca_rest))) / (1 + p.lambda_buf)         
    du[idx['nNOS']] = p.NOswitch_NE * p.V_maxNOS * a.CaM / (p.K_actNOS + a.CaM) - p.mu2_n * v.nNOS 
    du[idx['NO_n']] = a.p_NO_n - a.c_NO_n + a.d_NO_n          
    
    du[idx['K_p']] = a.J_BK_k / (p.VR_pa) + a.J_KIR_i / p.VR_ps - p.R_decay * (v.K_p - p.K_p_min) 
    du[idx['Ca_p']] = a.J_TRPV_k/ p.VR_pa + a.J_VOCC_i / p.VR_ps - p.Ca_decay_k * (v.Ca_p - p.Capmin_k)
    du[idx['v_k']] = p.gamma_i * ( -a.J_BK_k - a.J_K_k - a.J_Cl_k - a.J_NBC_k - a.J_Na_k - a.J_NaK_k - 2*a.J_TRPV_k - a.J_GABA_k)
    du[idx['Na_k']] = -a.J_Na_k - 3*a.J_NaK_k + a.J_NKCC1_k + a.J_NBC_k
    du[idx['K_k']] = -a.J_K_k + 2*a.J_NaK_k + a.J_NKCC1_k + a.J_KCC1_k - a.J_BK_k
    du[idx['HCO3_k']] = 2*a.J_NBC_k
    du[idx['Ca_k']] = a.B_cyt * (a.J_IP3 - a.J_pump + a.J_ER_leak - a.J_TRPV_k/p.r_buff)  
    du[idx['Cl_k']] = du[idx['Na_k']] + du[idx['K_k']] - du[idx['HCO3_k']] + p.z_Ca * du[idx['Ca_k']]
    
    du[idx['K_s']] = 1/p.VR_sa * (a.J_K_k - 2 * a.J_NaK_k - a.J_NKCC1_k - a.J_KCC1_k) + a.J_K_NEtoSC
    du[idx['Na_s']] = 1/p.VR_sa * (a.J_Na_k + 3 * a.J_NaK_k - a.J_NKCC1_k - a.J_NBC_k) + a.J_Na_NEtoSC
    du[idx['HCO3_s']] = 1/p.VR_sa * (-2*a.J_NBC_k)
    du[idx['w_k']] = a.phi_w * (a.w_inf - v.w_k)
    du[idx['I_k']] =  p.r_h * a.G - p.k_deg * v.I_k
  
    du[idx['h_k']] = p.k_on * (p.K_inh - (v.Ca_k + p.K_inh) * v.h_k)
    du[idx['s_k']] = -(a.B_cyt * (a.J_IP3 - a.J_pump + a.J_ER_leak)) / (p.VR_ER_cyt)
    du[idx['m_k']] =  p.trpv_switch * ((a.minf_k - v.m_k) / p.t_TRPV_k)
    du[idx['eet_k']] = p.V_eet * max(v.Ca_k - p.Ca_k_min, 0) - p.k_eet * v.eet_k
    du[idx['NO_k']] = a.p_NO_k - a.c_NO_k + a.d_NO_k
    du[idx['AA_k']] = (p.AA_m * p.AA_max)/(p.AA_m + max(v.Ca_k - p.Ca0,0))**2 * du[idx['Ca_k']] + (v.AA_i - v.AA_k)/p.tau_AA 
    
    du[idx['Ca_i']] = a.J_IP3_i - a.J_SR_uptake_i - a.J_extrusion_i + a.J_SR_leak_i - a.J_VOCC_i + a.J_CICR_i + a.J_NaCa_i - 0.1*a.J_stretch_i + a.J_Ca_coup_i
    du[idx['s_i']] = a.J_SR_uptake_i - a.J_CICR_i - a.J_SR_leak_i
    du[idx['v_i']] = p.gamma_i * ( -a.J_NaK_i - a.J_Cl_i - 2*a.J_VOCC_i - a.J_NaCa_i - a.J_K_i - a.J_stretch_i - a.J_KIR_i - a.J_GABA_i) + a.V_coup_i
    du[idx['w_i']] = p.lambda_i * (a.K_act_i - v.w_i) 
    du[idx['I_i']] = a.J_IP3_coup_i - a.J_degrad_i
    du[idx['NO_i']] = a.p_NO_i - a.c_NO_i + a.d_NO_i
    du[idx['E_b']] = -p.k1 * v.E_b * v.NO_i + p.k_1 * v.E_6c + a.k4 * a.E_5c     
    du[idx['E_6c']] = p.k1 * v.E_b * v.NO_i - (p.k_1 + p.k2) * v.E_6c - p.k3 * v.E_6c * v.NO_i
    du[idx['cGMP_i']] = p.V_max_sGC * a.E_5c - a.V_max_pde * v.cGMP_i / (p.K_m_pde + v.cGMP_i) 
    du[idx['H_i']] = a.f_NO * p.V_a * v.AA_i /(p.K_a + v.AA_i) +  p.V_f * v.AA_i / (p.K_f + v.AA_i) - p.lambda_h * v.H_i 
    du[idx['AA_i']] = (v.AA_k - v.AA_i)/p.tau_AA 
    
    du[idx['Ca_j']] = a.J_IP3_j - a.J_ER_uptake_j + a.J_CICR_j - a.J_extrusion_j + a.J_ER_leak_j + a.J_cation_j + p.J_0_j - a.J_stretch_j - a.J_Ca_coup_i
    du[idx['s_j']] = a.J_ER_uptake_j - a.J_CICR_j - a.J_ER_leak_j 
    du[idx['v_j']] = -1/p.C_m_j * (a.J_K_j + a.J_R_j) - a.V_coup_i
    du[idx['I_j']] = p.J_PLC - a.J_degrad_j - a.J_IP3_coup_i
    du[idx['eNOS']] = p.gam_eNOS * a.Act_eNOS_Ca * p.NOswitch_EC_CA  + (1 - p.gam_eNOS) * a.Act_eNOS_wss * p.NOswitch_EC_WSS - p.mu2_j * v.eNOS 
    du[idx['NO_j']] = a.p_NO_j - a.c_NO_j + a.d_NO_j 
    
    du[idx['Mp']] = p.wallMech * ( p.K_4 * v.AMp + a.K_1 * a.M - (a.K_2 + p.K_3) * v.Mp )
    du[idx['AMp']] = p.wallMech  * ( p.K_3 * v.Mp + a.K_6 * v.AM - (p.K_4 + a.K_5) * v.AMp )
    du[idx['AM']] = p.wallMech  * ( a.K_5 * v.AMp - (p.K_7 + a.K_6) * v.AM )
    du[idx['R']] =  p.R_init / p.eta_R * ( v.R * p.trans_p / a.h - a.E * (v.R - a.R_0) / a.R_0)
    
#    du[idx['GABA']] = - a.kappa_GABA * v.GABA + a.Ginput
    du[idx['GABA']] = - a.kappa_GABA * (v.GABA - p.GABAbase ) + p.alpha_GABA * (v.I_t - p.Imin) / (p.I_relative - p.Imin)
    du[idx['NPY']] = - p.beta_NPY * (v.NPY - p.NPYbase) + p.alpha_NPY * (v.I_t - p.Imin) / (p.I_relative - p.Imin)
    du[idx['Glu']] = - p.beta_Glu * v.Glu + p.GluSwitch * (a.f_Ke + a.f_GABA) 
        
    # Print to console the current time and percentage completed
    time_elapsed = time.time() - start_time
    print('\r-- Time {0:.2f} / {1:.0f} -- {2:.1f}% complete -- {3:.0f} sec elapsed --'.format( t/1e3, p.Tend/1e3, t/p.Tend*100, time_elapsed ), end='')  
    
    return du