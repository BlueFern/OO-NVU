# -*- coding: utf-8 -*-
"""
Algebraic variables
"""

import numpy as np
from state_variable_names import initialise_var_names
import model_functions as f
import parameters as p

# Define a class that will contain attributes = algebraic variables
class alg_vars():   
    pass

# Returns a class 'a' containing all algebraic variables
# u: 2D array of state variables vs time
# t: time vector
# idx: dictionary containing the state variable names ('keys') and their corresponding index ('values')
# input_data: experimental data used for the neuronal input functions
# total_waveform: an array containing the fast triangular input pulses when using the thalamic input T(t)
# location: for use inside a function ('in function') or an array for use outside of a function ('out function')
def set_algebraic_variables(u, t, idx, input_data, total_waveform, location):
    
    # Initialise variable names for clarity 
    v = initialise_var_names(u, idx, location)
    
    # Initialise class to contain the variables
    a = alg_vars()
       
    #%%
    '''Neuron'''

    if p.InputCase == 'ThalamicTrianglesZheng' or p.InputCase == 'ThalamicSquarePulse' or p.InputCase == 'ThalamicTriangles': # Thalamic input, Pinto et al
        
        a.T = f.input_T(t, input_data, total_waveform, location)
        a.arg_e = p.wee * v.E_t - p.wie * v.I_t + p.wte * a.T 
        a.arg_i = p.wei * v.E_t - p.wii * v.I_t + p.wti * a.T 
        a.f_arg_e = 5.12 / ( 1 + np.exp(-(a.arg_e - 15) / 4.16) )
        a.f_arg_i = 11.61 / ( 1 + np.exp(-(a.arg_i - 15) / 3.94) )
      
    else: # Normal P,Q input
        
        a.k_e = 1 / (1 + np.exp(-(p.a_e * (p.k_e_input - p.theta_e)))) - 1 / (1 + np.exp(p.a_e * p.theta_e))
        a.k_i = 1 / (1 + np.exp(-(p.a_i * (p.k_i_input - p.theta_i)))) - 1 / (1 + np.exp(p.a_i * p.theta_i))
    
        a.P = f.input_P(t, input_data, location)
        a.Q = f.input_Q(t, input_data, location)
    
        a.arg_e = p.c1 * v.E_t - p.c2 * v.I_t + a.P  # [-]           Excitatory cells response function input
        a.arg_i = p.c3 * v.E_t - p.c4 * v.I_t + a.Q  # [-]           Inhibitory cells response function input
    
        a.s_arg_e = 1 / (1 + np.exp(-(p.a_e * (a.arg_e - p.theta_e)))) - 1 / (1 + np.exp(p.a_e * p.theta_e))
        a.s_arg_i = 1 / (1 + np.exp(-(p.a_i * (a.arg_i - p.theta_i)))) - 1 / (1 + np.exp(p.a_i * p.theta_i))

    a.O2_p            = p.O2_0 * (1 - p.O2switch) + v.O2 * p.O2switch  # [mM]  Oxygen concentration dependent on ATP                                   
    a.J_pump2         = 2 * (1 + p.O2_0 / (((1 - p.alpha_O2) * a.O2_p)  + p.alpha_O2 * p.O2_0))**(-1) # [mM s**-1]       
    a.J_pump2_0       = 0.0952           # [mM s**-1]                      
    a.J_pump2_O2_0    = 1                # [mM s**-1]                         
    a.P_02            = (a.J_pump2 - a.J_pump2_0 ) / ( a.J_pump2_O2_0 - a.J_pump2_0)   # [-]  Normalised pump rate of oxygen            
    
    a.CBF             = p.CBF_init * (v.R**4 / p.R_init**4)             # [-] Normalised cerebral blood flow
                                
    a.J_O2_background = p.J_0 * a.P_02 * (1 - p.gamma_O2)  # [mM s**-1]    Background oxygen consumption                                                                                 
    a.J_pump1_sa      = (1 + (p.K_init_e / v.K_e))**(-2) * (1 +  (p.Na_init_sa / v.Na_sa)) ** (-3) # [mM s**-1]  Oxygen consumption in soma/axon pump
    a.J_pump1init_sa  = 0.0312           # [mM s**-1]   Initial consumption in soma/axon pump    
    a.J_pump1_d       = (1 + (p.K_init_e / v.K_e))**(-2) * (1 + (p.Na_init_d / v.Na_d))**(-3)   # [mM s**-1]  Oxygen consumption in dendrite pump
    a.J_pump1init_d   = 0.0312           # [mM s**-1]   Initial consumption in dendrite pump             
    a.J_O2_pump       = p.J_0 *  a.P_02 * p.gamma_O2 * (( a.J_pump1_sa +  a.J_pump1_d) / ( a.J_pump1init_sa +  a.J_pump1init_d)) # [mM s**-1]    Total oxygen consumption of the NaK pumps  
    
    # BOLD
#    a.f_out           = v.CBV**(1/p.d) + p.tau_TAT * (1/(p.tau_MTT + p.tau_TAT) * (  a.CBF/p.CBF_init  - v.CBV**(1/p.d) ))  # [-]   Buxton balloon model outflow                                     
#                                                   
    a.CaM = v.Ca_n / p.m_c                 # [uM]        
    a.p_NO_n = ( v.nNOS * p.V_max_NO_n * p.O2_n /  (p.K_mO2_n + p.O2_n) * p.LArg_n / (p.K_mArg_n + p.LArg_n) )
    a.c_NO_n = p.k_O2_n * v.NO_n**2 * p.O2_n   
    a.tau_nk = p.x_nk ** 2 /  (2 * p.D_cNO)                          # [ms]
    a.d_NO_n = (v.NO_k - v.NO_n) /  a.tau_nk         
    
    # Flux from ECS to SC
    a.dKedt = -p.beta_K_e * (v.K_e - p.K_eBase) + (p.alpha_K_e * p.beta_K_e) * ((abs(v.E_t - v.I_t) - p.EImin)/ (p.EI_relative - p.EImin))          # K_e derivative as neuron activation in Ostby            
    
    a.J_K_NEtoSC = p.k_syn * a.dKedt * 1e3                                   # Input flux: K+ flux into SC
    a.J_Na_NEtoSC = -p.k_syn * a.dKedt * 1e3                                 # Input flux: Na+ flux out of SC
    
    #%%
    ''' Astrocyte'''
    
    # Electroneutrality condition
    a.Cl_s = v.Na_s + v.K_s - v.HCO3_s
    
    # Nernst potentials (in mV)
    a.E_K_k = p.ph / p.z_K * np.log(v.K_s / v.K_k)
    a.E_Na_k = p.ph / p.z_Na * np.log(v.Na_s / v.Na_k)
    a.E_Cl_k = p.ph / p.z_Cl * np.log(a.Cl_s / v.Cl_k)
    a.E_NBC_k = p.ph / p.z_NBC * np.log((v.Na_s * v.HCO3_s**2) / (v.Na_k * v.HCO3_k**2))
    a.E_BK_k = p.ph / p.z_K * np.log(v.K_p / v.K_k)        
    a.E_TRPV_k = p.ph / p.z_Ca * np.log(v.Ca_p / v.Ca_k)                   # Nernst potential TRPV

    # Flux through the Sodium Potassium pump
    a.J_NaK_k = p.J_NaK_max * v.Na_k**1.5 / (v.Na_k**1.5 + p.K_Na_k**1.5) * v.K_s / (v.K_s + p.K_K_s)
    
    # Fluxes
    a.J_BK_k = p.G_BK_k * v.w_k * (v.v_k - a.E_BK_k)                      # BK flux (uM/ms)
    a.J_BK_p = a.J_BK_k / p.VR_pa                                     # K+ influx into the PVS (uM/ms)
    a.J_K_k = p.G_K_k * (v.v_k - a.E_K_k)
    a.J_Na_k = p.G_Na_k * (v.v_k - a.E_Na_k)
    a.J_NBC_k = p.G_NBC_k * (v.v_k - a.E_NBC_k)
    a.J_KCC1_k = p.G_KCC1_k * p.ph * np.log((v.K_s * a.Cl_s) / (v.K_k * v.Cl_k))
    a.J_NKCC1_k = p.G_NKCC1_k * p.ph * np.log((v.Na_s * v.K_s * a.Cl_s**2) / (v.Na_k * v.K_k * v.Cl_k**2))
    a.J_Cl_k = p.G_Cl_k * (v.v_k - a.E_Cl_k)
    
    # Calcium Equations
    # Flux
    a.J_IP3 = p.J_max * ( v.I_k / (v.I_k + p.K_I) *  v.Ca_k / (v.Ca_k + p.K_act) * v.h_k)**3 * (1 - v.Ca_k / v.s_k)
    a.J_ER_leak = p.P_L * (1 - v.Ca_k / v.s_k)
    a.J_pump = p.V_max * v.Ca_k**2 / (v.Ca_k**2 + p.k_pump**2)
    a.J_TRPV_k = p.G_TRPV_k * v.m_k * (v.v_k - a.E_TRPV_k)
    
    
    # Other equations
    a.B_cyt = 1 / (1 + p.BK_end + p.K_ex * p.B_ex / (p.K_ex + v.Ca_k)**2)
    
    a.v_3 = p.v_6 - p.v_5 / 2 * np.tanh((v.Ca_k - p.Ca_3) / p.Ca_4)
    
    # Parent Calcium equations 
    a.w_inf = 0.5 * (1 + np.tanh((v.v_k + p.eet_shift * v.eet_k - a.v_3) / p.v_4))
    a.phi_w = p.psi_w * np.cosh((v.v_k - a.v_3) / (2 * p.v_4))
    
    # TRPV Channel open probability equations
    a.H_Ca_k = v.Ca_k / p.gam_cai_k + v.Ca_p / p.gam_cae_k
    a.eta = (v.R - p.R_init) / (p.R_init)
    a.minf_k = (1 / (1 + np.exp(-(a.eta - p.epshalf_k) / p.kappa_k))) * ((1 / (1 + a.H_Ca_k)) * (a.H_Ca_k + np.tanh((v.v_k - p.v1_TRPV_k) / p.v2_TRPV_k))) 
    
    # NO pathway
    a.tau_ki = p.x_ki ** 2 /  (2 * p.D_cNO)                          # [ms]
    a.p_NO_k = 0
    a.c_NO_k = p.k_O2_k * v.NO_k**2 * p.O2_k                           # [uM/ms]
    a.d_NO_k = (v.NO_n - v.NO_k) / a.tau_nk + (v.NO_i - v.NO_k) / a.tau_ki     # [uM/ms]

    #%%    
    ''' SMC/EC'''
    
    # SMC fluxes
    a.J_IP3_i = p.F_i * v.I_i**2 / (p.K_r_i**2 + v.I_i**2)           # IP3 channel
    a.J_SR_uptake_i = p.B_i * v.Ca_i**2 / (p.c_b_i**2 + v.Ca_i**2)   # SERCA pump
    a.J_CICR_i = p.C_i * v.s_i**2 / (p.s_c_i**2 + v.s_i**2) * v.Ca_i**4 / (p.c_c_i**4 + v.Ca_i**4)
    a.J_extrusion_i = p.D_i * v.Ca_i * (1 + (v.v_i - p.v_d) / p.R_d_i)
    a.J_SR_leak_i = p.L_i * v.s_i
    a.J_NaCa_i = p.G_NaCa_i * v.Ca_i / (v.Ca_i + p.c_NaCa_i) * (v.v_i - p.v_NaCa_i)
    
    a.h = 0.1 * v.R 
    a.J_stretch_i = p.G_stretch / (1 + np.exp(-p.alpha_stretch*(p.trans_p_mmHg*v.R/a.h - p.sigma_0))) * (v.v_i - p.E_SAC)
    
    a.J_Cl_i = p.G_Cl_i * (v.v_i - p.v_Cl_i)
    a.J_NaK_i = p.F_NaK_i
    a.J_K_i   = p.G_K_i * v.w_i * (v.v_i - p.v_K_i)
    a.J_degrad_i = p.k_d_i * v.I_i
    
    a.v_KIR_i = p.z_1 * v.K_p - p.z_2
    a.G_KIR_i = p.F_KIR_i * np.exp(p.z_5 * v.v_i + p.z_3 * v.K_p)           # Changed by including np.exp(-z_4) parameter into F_KIR_i parameter, no large constant in np.exponential equation
    a.J_KIR_i = a.G_KIR_i * (v.v_i - a.v_KIR_i)
    
    # EC fluxes
    a.J_IP3_j = p.F_j * v.I_j**2 / (p.K_r_j**2 + v.I_j**2)
    a.J_ER_uptake_j = p.B_j * v.Ca_j**2 / (p.c_b_j**2 + v.Ca_j**2)  
    a.J_CICR_j = p.C_j * v.s_j**2 / (p.s_c_j**2 + v.s_j**2) * v.Ca_j**4 / (p.c_c_j**4 + v.Ca_j**4)
    a.J_extrusion_j = p.D_j * v.Ca_j
    a.J_stretch_j = p.G_stretch / (1 + np.exp(-p.alpha_stretch*(p.trans_p_mmHg*v.R/a.h - p.sigma_0))) * (v.v_j - p.E_SAC)
    a.J_ER_leak_j = p.L_j * v.s_j
    a.J_cation_j = p.G_cat_j * (p.E_Ca_j - v.v_j) * 0.5 * (1 + np.tanh((np.log10(v.Ca_j) - p.m_3_cat_j) / p.m_4_cat_j))
    a.J_BK_Ca_j = 0.2 * (1 + np.tanh( ((np.log10(v.Ca_j) - p.c) * (v.v_j - p.bb_j) - p.a_1_j) / (p.m_3b_j * (v.v_j + p.a_2_j*(np.log10(v.Ca_j) - p.c) - p.bb_j)**2 + p.m_4b_j)))
    a.J_SK_Ca_j = 0.3 * (1 + np.tanh((np.log10(v.Ca_j) - p.m_3s_j) / p.m_4s_j))
    a.J_K_j = p.G_tot_j * (v.v_j - p.v_K_j) * (a.J_BK_Ca_j + a.J_SK_Ca_j)
    a.J_R_j = p.G_R_j * (v.v_j - p.v_rest_j)
    a.J_degrad_j = p.k_d_j * v.I_j
    
    # Coupling
    a.V_coup_i = -p.G_coup * (v.v_i - v.v_j)
    a.J_IP3_coup_i = -p.P_IP3 * (v.I_i - v.I_j)
    a.J_Ca_coup_i = -p.P_Ca * (v.Ca_i - v.Ca_j)
    
    a.c_w_i = 1/2 * (1 + np.tanh((v.cGMP_i - 10.75)/0.668) )
    a.K_act_i = (v.Ca_i + a.c_w_i)**2 / ((v.Ca_i + a.c_w_i)**2 + p.alpha_act_i * np.exp(-(v.v_i - p.v_Ca3_i - p.Hshift*(v.H_i-p.H0)) / p.R_K_i)) # Is shifted by 20-HETE ******** 
    a.tau_wss = v.R/2 * p.delta_p_L # from Dormanns 2016
    
    # NO pathway 
    a.tau_ki = p.x_ki ** 2 /  (2 * p.D_cNO)                          # [ms]
    a.tau_ij = p.x_ij ** 2 /  (2 * p.D_cNO)                          # [ms]
    a.p_NO_i = 0
    a.c_NO_i = p.k_dno * v.NO_i
    a.d_NO_i = (v.NO_k - v.NO_i) / a.tau_ki + (v.NO_j - v.NO_i) / a.tau_ij     # [uM]
    a.k4 = p.C_4 * v.cGMP_i**2
    a.E_5c = 1 - v.E_b - v.E_6c
    a.V_max_pde = p.k_pde * v.cGMP_i
    a.R_cGMP2 = v.cGMP_i**2 / (v.cGMP_i**2 + p.K_m_mlcp**2)
    
    a.O2_j = v.O2*1e3  # Oxygen in EC taken as O2 from lumen (diffusion very fast so plausible!) instead of constant, in uM
    a.p_NO_j = ( p.V_NOj_max * v.eNOS * a.O2_j / (p.K_mO2_j + a.O2_j) * p.LArg_j / (p.K_mArg_j + p.LArg_j) )
    a.c_NO_j = p.k_O2 * v.NO_j**2 * a.O2_j# consumption by oxygen 
    
    a.J_lumen = - v.NO_j * 4 * p.D_cNO / (25**2) 
    a.d_NO_j = (v.NO_i - v.NO_j) / a.tau_ij + a.J_lumen 
    
    a.W_wss = p.W_0 * (a.tau_wss + np.sqrt(16 * p.delta_wss**2 + a.tau_wss**2) - 4 * p.delta_wss)**2 / (a.tau_wss + np.sqrt(16 * p.delta_wss**2 + a.tau_wss**2))    
    a.F_wss = 1 / (1 + p.alp * np.exp(- a.W_wss)) - 1 / (1 + p.alp) 
    a.Act_eNOS_Ca = p.K_dis * v.Ca_j / (p.K_eNOS + v.Ca_j) 
    a.Act_eNOS_wss = p.g_max * a.F_wss  
    
    # 20-HETE
    a.f_NO = 1 / (1 + np.exp((v.NO_i - p.NO_rest) / p.R_NO))
    
    #%%    
    ''' Wall mechanics'''
    
    a.K_1 = p.gamma_cross * v.Ca_i**p.n_cross                          
    a.K_6 = a.K_1                                                      
    a.K_2 = 58.1395 * p.k_mlcp_b + 58.1395 * p.k_mlcp_c * a.R_cGMP2    
    a.K_5 = a.K_2                                                      
    
    a.M = 1 - v.AM - v.AMp - v.Mp        
    
    # Mechanical Equations
    a.F_r = v.AMp + v.AM
    a.E = p.E_passive + a.F_r * (p.E_active - p.E_passive)
    a.R_0 = p.R_init + a.F_r * (p.alpha - 1) * p.R_init
    
    #%%    
    ''' Oxygen extraction fraction & CMRO2 equations'''
     
    a.f_in_dim = a.CBF/(p.CBF_init)                                           # inflow of blood (normalised CBF)
    a.f_in = a.f_in_dim/p.f_in0                                                           
    a.f_out = v.CBV**(1/p.d) + p.tau_TAT * ( 1/(p.tau_MTT + p.tau_TAT) * ( a.f_in  - v.CBV**(1/p.d) ) )  # [-]   Buxton balloon model outflow                                     
    a.OEF = 1 - (1 - p.E_0)**(1/a.f_in)                                       # oxygen extraction fraction based on Buxton 1997    
    a.CMRO2 = a.f_in * a.OEF / p.E_0                                               # Cerebral metabolic rate of oxygen consumption
    a.J_O2_vascular   = a.CBF * a.OEF / p.E_0                                      # [mM s**-1] Vascular supply of oxygen, previously CBF * ((p.O2_b - O2) / (p.O2_b - p.O2_0))   
 
    #%%    
    ''' GABA stuff'''   
    # GABA-T activity, = 1 with normal levels of NO and = 2 when NO=0 (NO inhibits GABA-T activity) [-]
#    a.G_Tact = p.G_Tmin + (p.G_Tmax - p.G_Tmin) * np.exp(-p.p2 * v.NO_n) # old exponential form
    a.G_Tact = 0.5 * ( (p.G_Tmax + p.G_Tmin) - (p.G_Tmax - p.G_Tmin) * np.tanh( (v.NO_n - p.GT_midpoint)/p.GT_slope ) ) # new sigmoid form
    
    a.kappa_GABA = p.beta_GABA * a.G_Tact                                     # degradation of GABA due to GABA-T [s**-1]
    a.g_GABA = p.G_GABA * 0.5 * ( 1 + np.tanh( (v.GABA - p.g_midpoint) / p.g_slope ) ) # Conductance of GABA dependent Cl channel, where conductance=0 when GABA=0 and conductance=G_Cl (value taken from SMC model) when GABA is max
    
    # Fluxes through the GABA activated Cl channels
    a.J_GABA_k = a.g_GABA * (v.v_k - p.E_GABA)
    a.J_GABA_i = a.g_GABA * (v.v_i - p.E_GABA)
    
    # Glutamate algebraic variables
    a.f_Ke = p.beta_Glu * 0.5 * ( 1 + np.tanh( (v.K_e - p.Ke_switch) / p.Glu_slope) )  # Release of Glu from the excitatory neuron during stimulation (occurs when Ke > 5 mM) [-]
    a.f_GABA = a.kappa_GABA * v.GABA       # Production of Glu from the degradation of GABA (GABA -> succinic semialdehyde + Glu) [s**-1]
    
    a.w_NR2A = v.Glu / (p.K_mA + v.Glu)     # [-] 
    a.w_NR2B = v.Glu / (p.K_mB + v.Glu)     # [-] 
    a.I_Ca = (-4 * p.v_n * p.G_M * p.P_Ca_P_M * (p.Ca_ex / p.M)) / (1 + np.exp(-0.08 * (p.v_n + 20))) * (np.exp(2 * p.v_n / p.ph)) / (1 - np.exp(2 * p.v_n / p.ph))   # [fA]                                    
    a.I_Ca_tot = a.I_Ca * (p.n_NR2A * a.w_NR2A + p.n_NR2B * a.w_NR2B)   # [fA]  
    
    a.rho = p.rho_min + (p.rho_max - p.rho_min) * v.Glu
    a.G = (a.rho + p.delta) / (p.K_G + a.rho + p.delta)

    #%%    
    ''' NPY: released in interneurons, has a vasoconstrictive effect by opening the VOCCs'''
    # Conductance increases from baseline when NPY increases
    a.g_VOCC = p.G_Ca_i * (1 + p.NPYswitch * p.npy_increase * 0.5 * ( 1 + np.tanh( (v.NPY - p.npy_midpoint) / p.npy_slope ) ) )
    a.J_VOCC_i = a.g_VOCC * (v.v_i - p.v_Ca1_i) / (1 + np.exp(-(v.v_i - p.v_Ca2_i) / p.R_Ca_i)) 
    
    return a



