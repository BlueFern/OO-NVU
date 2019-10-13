# -*- coding: utf-8 -*-
"""
Model parameters

This file can be loaded as a module (e.g. put "import parameters as p" then access parameters as p.startpulse, p.dt etc)
"""

# Latest changes:   tau_MTT 3e3 -> 1e3
#                   k_syn 10.5 -> 11.5
#                   G_GABA .2 * G_Cl_i -> .24 * G_Cl_i
#                   npy_increase 0.04 -> 0.06
#                   Hshift 30 -> 10
# Changed to fit with correctly scaled whisker/opto data (previously HbO and HbR were divided by 100, now HbO is divided by 70 and HbR by 30)

''' Different model cases'''
InputCase = 'SquarePulse'            # 'SquarePulse': Model with P,Q as a square input pulse [default]
                                     # 'ThalamicSquarePulse': Thalamic input T (Pinto et al) as a square input pulse
                                     # 'Inhibitory': model with inhibitory neurons only (P = 0, Q as a square pulse)
                                     # 'ZhengData': Model with P,Q scaled by Zheng experimental data
                                     # '2SquarePulses': Model with 2 square input pulses following the Zheng paradigm (second pulse 2 sec long)
                                     # 'ZhengFittedParams': Model with P,Q scaled by Zheng experimental data and several parameters fitted to experiment (Yves)
                                     # 'ThalamicTriangles': Thalamic input T (Pinto et al) as a series of fast triangular pulses determined in T_pulses() (in model_functions module)
                                     # 'ThalamicTrianglesZheng': Thalamic input T (Pinto et al) as a series of fast triangular pulses determined in T_pulses() (in model_functions module) scaled by Zheng experimental data

NeuronType = 'whisker'   # 'whisker': normal excitatory neuron (whisker stimulation, somatosensory cortex) [default]
                         # 'opto': inhibitory interneuron (optogenetic response), note that this isn't working properly for the Thalamic input cases yet!!
                         # 'whiskerandopto': combo of both whisker and opto

NOswitch = 'normal'     # 'normal': eNOS and nNOS production [default]
                        # '7NI': no nNOS production
                        # 'LNAME': no eNOS or nNOS production
                        
if (NeuronType == 'whisker' or NeuronType == 'whiskerandopto') and InputCase != 'Inhibitory':
    Pinput = 1.0 
    Qinput = 1.0
else: 
    Pinput = 0.0 # inhibitory input only, for NeuronType = 'opto' or special 'Inhibitory' whisker case
    Qinput = 1.0
                    
if NOswitch == 'normal':
    NOswitch_NE = 1      # Turn on neuronal NO production 
    NOswitch_EC_WSS = 1  # Turn on WSS induced eNOS production 
    NOswitch_EC_CA = 1   # Turn on Ca2+ induced eNOS production 
elif NOswitch == '7NI':
    NOswitch_NE = 0      # Turn off neuronal NO production 
    NOswitch_EC_WSS = 1  # Turn on WSS induced eNOS production 
    NOswitch_EC_CA = 1   # Turn on Ca2+ induced eNOS production 
elif NOswitch == 'LNAME':
    NOswitch_NE = 0      # Turn off neuronal NO production 
    NOswitch_EC_WSS = 0  # Turn off WSS induced eNOS production 
    NOswitch_EC_CA = 0   # Turn off Ca2+ induced eNOS production

''' Input [ms] '''
startpulse = 500e3       # Start time of input stimulation, set as 500 sec to allow the model to reach a steady state
lengthpulse = 2e3        # Length of stimulation
Tend = 550e3             # End of simulation
dt = 0.001e3             # Time step (this value is overwritten to be small if InputCase = 'ThalamicTriangles' or 'ThalamicTrianglesZheng') [0.001e3 default]
double_pulse = 0         # 1: use the second pulse in the Zheng input data, 0: only one pulse. Used for InputCase = 'ZhengData' or 'ZhengFittedParams' [0:default]

''' Neuron '''

if NeuronType == 'whisker': # less GABA and NPY as there is a smaller population of inhibitory neurons in the whisker case
    alpha_GABA = 2.6e-3     # used to scale inflow of GABA due to increase in I_t, M.E. (Maite)
    alpha_NPY = 2.6e-3      # used to scale inflow of NPY due to increase in I_t, M.E. (Maite)
elif NeuronType == 'opto' or NeuronType == 'whiskerandopto':
    alpha_GABA = 4.2e-3     # used to scale inflow of GABA due to increase in I_t, M.E. (Maite)
    alpha_NPY = 4.2e-3      # used to scale inflow of NPY due to increase in I_t, M.E. (Maite)    

if InputCase == 'ThalamicSquarePulse' or InputCase == 'ThalamicTriangles' or InputCase == 'ThalamicTrianglesZheng': # Thalamic input (Pinto et al)
    EI_relative = 1.499         # [-]    Maximum value for abs(E-I) after external inputs (block pulse)
    EImin = 0.139               # [-]    Minimum value for abs(E-I)
    I_relative = 1.512          # [-]    Maximum value for I after external inputs (block pulse)
    Imin = 0.225                # [-]    Minimum value for I
else: # P, Q input
    EI_relative = 0.268         # [-]    Maximum value for abs(E-I) after external inputs
    EImin = 0                   # [-]    Minimum value for abs(E-I) 
    Imin = 0                    # [-]    Minimum value for I    
    if NeuronType == 'whisker' or NeuronType == 'whiskerandopto':
        I_relative = 0.179      # [-]    Maximum value for I after external inputs 
    elif NeuronType == 'opto':
        I_relative = 0.016      # [-]    Maximum value for I after external inputs, smaller maximum when only inhibitory neurons active (P=0)

# Wilson and Cowan model parameters for K_e, Na_sa, Na_d
# Alpha, beta and rest values adjusted for the different inputs
if InputCase == 'SquarePulse' or InputCase == '2SquarePulses' or InputCase == 'Inhibitory':
    alpha_K_e = 2.0
    beta_K_e = 4.2e-3
    alpha_Na_sa = 4.23
    beta_Na_sa = 0.39e-3
    alpha_Na_d = -2.12
    beta_Na_d = 0.75e-3
    K_eBase = 3.5
    Na_saBase = 9.37
    Na_dBase = 9.42
    tau_e = 3.0      
    tau_i = 3.0           
elif InputCase == 'ZhengData':
    alpha_K_e = 2.2
    beta_K_e = 4.2e-3
    alpha_Na_sa = 4.23
    beta_Na_sa = 0.65e-3
    alpha_Na_d = -2.12
    beta_Na_d = 0.7e-3
    K_eBase = 3.567
    Na_saBase = 9.237
    Na_dBase = 9.318
    tau_e = 3.0        
    tau_i = 3.0    
elif InputCase == 'ThalamicSquarePulse' or InputCase == 'ThalamicTriangles' or InputCase == 'ThalamicTrianglesZheng': # Thalamic input
    alpha_K_e = 1.8
    beta_K_e = 4.2e-3 #1.7e-3
    alpha_Na_sa = 4.23
    beta_Na_sa = 0.39e-3
    alpha_Na_d = -2.12
    beta_Na_d = 0.75e-3
    K_eBase = 3.5
    Na_saBase = 9.37
    Na_dBase = 9.42
    tau_e = 5.0        
    tau_i = 15.0    
elif InputCase == 'ZhengFittedParams':
    alpha_K_e = 3.3    # 3.3 with wallMech = 2.0 for CBF fit to Zheng, 3.0 with wallMech = 2.0 for fit to Berwick
    beta_K_e = 10e-3
    alpha_Na_sa = 4.23
    beta_Na_sa = 0.65e-3
    alpha_Na_d = -2.12
    beta_Na_d = 0.7e-3
    K_eBase = 3.567
    Na_saBase = 9.237
    Na_dBase = 9.318
    tau_e = 3        
    tau_i = 3  

 # Max size of thalamic input function T(t)
if InputCase == 'ThalamicTrianglesZheng':   # multiple fast triangular pulses scaled by Zheng input data
    Tinput = 1
elif InputCase == 'ThalamicTriangles':      # multiple fast triangular pulses
    Tinput = 0.8
elif InputCase == 'ThalamicSquarePulse':    # square block pulse, effectively a time 'average' of a series of fast triangular pulses
    Tinput = 0.5
else:
    Tinput = 0
    
# Zheng 2010 paradigm
if InputCase == 'ZhengData' or InputCase == 'ThalamicTrianglesZheng' or InputCase == '2SquarePulses' or InputCase == 'ZhengFittedParams':
    ISI = 8e3                                           # Inter stimulus interval (ISI), i.e. the time in ms between stimulations [8e3:default]
    ISI_array = [0.6e3, 1e3, 2e3, 3e3, 4e3, 6e3, 8e3]   # array of possible ISI
    ISI_idx = ISI_array.index(ISI)                      # INDEX for the ISI
    
    # stim = INDEX for length of initial stimulation [2,8,16] beginning at zero
    if lengthpulse == 2e3:
        stim = 0
    elif lengthpulse == 8e3:
        stim = 1
    elif lengthpulse == 16e3:
        stim = 2
    else:
        print("When using this InputCase you must have lengthpulse = 2, 8 or 16 sec")

if InputCase == 'ZhengData' or InputCase == '2SquarePulses' or InputCase == 'ZhengFittedParams':
    secondpulse = startpulse + lengthpulse + ISI  
    secondlength = 2e3                       # length of second pulse
else:
    ISI = 0       # No second pulse
    secondlength = 0

if InputCase == 'ThalamicTriangles' or InputCase == 'ThalamicTrianglesZheng':
    dt = 0.001e3 # Thalamic input with fast triangles needs a small timestep

# E and I parameters
c1 = 12              # [-]
c2 = 10              # [-]
c3 = 13              # [-]
c4 = 11              # [-]
r_e = 1.0            # [ms]?
r_i = 4.0            # [ms]?
k_e_input = 10
k_i_input = 10
a_e = 1.2            # [-]
theta_e = 2.8        # [-]
a_i = 1.0            # [-]
theta_i = 4.0        # [-]

# Thalamic input parameters (Pinto et al)
wee = 42
wie = 24.6
wei = 42
wii = 18
wte = 53.43
wti = 68.4

# Glutamate parameters
Glu_slope = 0.1      # [mM] Slope of sigmoidal (M.E.)
Ke_switch = 5      # [mM] Threshold past which glutamate vesicle is released (M.E.)                               

# Elshin model constants
Farad = 96.485       # [C / mmol]
ph = 26.6995         # [RT/F]
k_syn = 11.5         # [-] Coupling scaling factor, values can range between 1.06 to 14.95 according to estimation from experimental data of Ventura and Harris . Original 11.5

# BOLD constants
d = 0.4              # [-]
a_1 = 3.4            # [-]   
a_2 = 1              # [-]
V_0 = 0.03           # [-]
tau_MTT = 1e3        # [ms]  3e3
tau_TAT = 5e3       # [ms] 

# Oxygen model
alpha_O2 =  0.05     # [-] Percentage of ATP production indepent of O2 0.05
K_init_e = 2.9       # [mM]
Na_init_sa = 10      # [mM]
Na_init_d = 10       # [mM]

CBF_init = 3.2e-2    # [-]3.2e-2 *****
gamma_O2 = 0.1       # [-] 0.1 

# NO pathway 
m_c = 4              # [-] 
K_mA = .352           # [uM] - fit to Santucci2008 650/1846
K_mB = 1.517          # [uM] - fit to Santucci2008 2800/1846
v_n = -40            # [mV]  the neuronal membrane potential , assumed to be approx constant in this model
G_M = 0.46           # original value 46e9 !! [pS] = [kg**-1*m**-2*ms**3*A**2] --> converted to ms
P_Ca_P_M = 3.6       # [-]  the relative conductance of the NMDA channel to Ca2+ compared to monovalent ions
Ca_ex = 2e3          # [uM]  the external calcium concentration (in Comerford+David2008: 1.5 mM!)
M = 1.3e5            # [uM]  the concentration of monovalent ions in the neuron
n_NR2A = 0.63        # [-]  average number of NR2A NMDA receptors per synapse (Santucci2008)
n_NR2B = 11          # [-]  average number of NR2B NMDA receptors per synapse (Santucci2008)
V_max_NO_n = 4.22e-3 # [ms**-1]  maximum catalytic rate of NO production (Chen2006) - obtained from fig 6 & equ 17 & 18
O2_n = 200           # [uM]  tissue O2 concentration in the neuron (M.E.)
K_mO2_n = 243        # [uM]  Chen2006
LArg_n = 100         # [uM]  
K_mArg_n = 1.5       # [uM]      
k_O2_n = 9.6e-9      # [uM**-2 ms**-1]  # (Kavdia2002)
x_nk = 25            # [um]   (M.E.)
D_cNO = 3300e-3      # [um**2 ms**-1]  Diffusion coefficient NO (Malinski1993)   
V_spine = 8e-5       # [pL]  volume of the neuronal dritic spine Santucci2008
k_ex = 1600e-3       # [ms**-1]  decay rate constant of internal calcium concentration Santucci2008
Ca_rest = 0.1        # [uM]  resting calcium concentration (in Comerford+David2008: 2.830 mM in Santucci2008P: 0.1 \muM)
lambda_buf = 20      # [-]  buffer capacity Santucci2008
V_maxNOS = 0.697e-3     # [uM ms**-1]  M.E. original values 25e-6
K_actNOS = 0.8   # [uM]  original values 9.27e-5
mu2_n = 0.2e-2    # [ms**-1]  original values 0.0167e-3 

# Switches
GluSwitch = 1        # Turn on/off glutamate input [1:default]
O2switch = 1         # O2switch = 0 ATP is plentiful, O2switch = 1 ATP is limited (oxygen-limited regime) [1:default]
LCpathway = 0        # Locus coeruleus pathway used for InputCase = 'ZhengData', 'ThalamicTrianglesZheng', 'ZhengFittedParams' [0:default]

''' Astrocyte '''
    
# Conductances in uM mV**-1 ms**-1
G_BK_k = 10.25e-3       
G_K_k = 6907.77e-3    
G_Na_k = 226.94e-3    
G_NBC_k = 130.74e-3   
G_KCC1_k = 1.728e-3
G_NKCC1_k = 9.568e-3 
G_Cl_k = 151.93e-3      

J_NaK_max = 2.3667e1                          # [uM ms**-1]    

# Switches to turn on and off some things
rhoSwitch = 1 
trpv_switch = 1 

rho_min = 0.1     
rho_max = 0.7  

v_4 = 8 #mV
v_5 = 15 #mV
v_6 = -55 #mV -55
Ca_3 = 0.4 # uM 
Ca_4 = 0.35 # uM 

# Scaling Constants
R_s = 2.79e-8 # m
R_k = 6e-8 # m
R_tot = 8.79e-8 # m
VR_sa = R_s/R_k

# Calcium in the Astrocyte Equations Constants
delta = 1.235e-2 
VR_ER_cyt = 0.185

J_max = 2880e-3                                  # [uM ms**-1] 
K_I = 0.03 # uM
K_act = 0.17 #uM
k_on = 2e-3                                      # [uM ms**-1]  
K_inh = 0.1 #uM
Hshift = 10  # previously 30 then 10

r_h = 4.8e-3                                     # [uM/ms] *****
k_deg = 1.25e-3                                  # [ms**-1]
V_eet = 72e-3                                    # [ms**-1]
k_eet = 7.2e-3                                   # [ms**-1]
Ca_k_min = 0.1 # uM
eet_shift = 2 #mV/uM


#TRPV4
Capmin_k = 2000 #uM
C_astr_k = 40#pF
gam_cae_k = 200 #uM
gam_cai_k = 0.01 #uM
epshalf_k = 0.1 # 
kappa_k = 0.1
v1_TRPV_k = 120 #mV
v2_TRPV_k = 13 #mV
t_TRPV_k = 0.9e3                                 # [ms] 
Ca_decay_k = 0.5e-3                              # [ms**-1]
G_TRPV_k = 3.15e-7                               # [ms**-1]
r_buff = 0.05 # Rate at which Ca2+ from the TRPV4 channel at the foot is buffered compared to rest of channels on the astrocyte body [-]

# Perivascular space
VR_pa = 0.001# [-]
VR_ps = 0.001# [-]
K_p_min = 3e3# uM
R_decay = 0.15e-3                           # [ms**-1]  0.15e-3   ***0.15e-3*3

# Fluxes Constants
R_g = 8.315 #J mol**-1 K**-1 Gas constant
T = 300 # K Temperature

K_Na_k = 10000 # uM
K_K_s = 1500 # uM

P_L = 0.0804e-3                                  # [uM ms**-1]
V_max = 20e-3                                    # [uM ms**-1]
k_pump = 0.24 #uM

# Additional Equations Astrocyte Constants
z_K = 1# [-]
z_Na = 1# [-]
z_Cl = -1# [-]
z_NBC = -1# [-]
z_Ca = 2# [-]
BK_end = 40# [-]
K_ex = 0.26 #uM
B_ex = 11.35 #uM
K_G = 8.82 
psi_w = 2.664e-3                                 # [ms**-1]

# NO Pathway
x_ki = 25            # [um]   (M.E.)
k_O2_k = 9.6e-9      # [uM**-2 ms**-1]   (Kavdia2002)
O2_k = 200           # [uM]   (M.E.)

''' SMC/EC '''

# Smooth Muscle Cell ODE Constants
gamma_i = 1970 #mV uM**-1
lambda_i = 45e-3                                 # [ms**-1]

# Endothelial Cell ODE Constants
C_m_j = 25.8                                     # pF veranderen naar iets met ms??
J_PLC = 0.11e-3                                 #0.11 for steady state, 0.3 for oscillations # [uM ms**-1] 
J_0_j = 0.029e-3                                 # [uM ms**-1] constant Ca influx (EC)

# Smooth Muscle Cell Flux Constants
F_i = 0.23e-3                                    # [uM ms**-1]      # IP3/RYR channel strength
K_r_i = 1 #uM

B_i = 2.025e-3                                   # [uM ms**-1]     # SERCA pump strength     
c_b_i = 1.0 #uM

C_i = 55e-3                                      # [uM ms**-1]
s_c_i = 2.0 #uM
c_c_i = 0.9 #uM

D_i = 0.24e-3                                    # [ms**-1]
v_d = -100 #mV
R_d_i = 250 #mV

L_i = 0.025e-3                                   # [ms**-1]

G_Ca_i = 1.29e-6                                 # [uM mV**-1 ms**-1]
v_Ca1_i = 100 #mV
v_Ca2_i = -24 #mV
R_Ca_i = 8.5 #mV

G_NaCa_i = 3.16e-6                               # [uM mV**-1 ms**-1]
c_NaCa_i = 0.5 #uM
v_NaCa_i = -30 #mV

G_stretch = 6.1e-6                               # [uM mV**-1 ms**-1]   (Also EC parameter)
alpha_stretch = 7.4e-3 # mmHg**-1     (Also EC parameter)
trans_p_mmHg = 30 # mmHg             (Also EC parameter) transmural pressure. 30 mmHg = 4000 Pa
sigma_0 = 500 # mmHg                 (Also EC parameter)
E_SAC = -18 # mV                     (Also EC parameter)

F_NaK_i = 4.32e-5                                # [uM ms**-1]

G_Cl_i = 1.34e-6                                 # [uM mV**-1 ms**-1]
v_Cl_i = -25 #mV

G_K_i = 4.46e-6                                  # [uM mV**-1 ms**-1]
v_K_i = -94 #mV

F_KIR_i = 1.285e-9#0.381                       # [uM mV**-1 ms**-1] F_KIR_i = 1.285e-6
k_d_i = 0.1e-3                                   # [ms**-1]

# othelial Cell Flux Constants
F_j = 0.23e-3                                    # [uM ms**-1]
K_r_j = 1 #uM
B_j = 0.5e-3                                     # [uM ms**-1]
c_b_j = 1 #uM
C_j = 5e-3                                       # [uM ms**-1]
s_c_j = 2 #uM
c_c_j = 0.9 #uM
D_j = 0.24e-3                                    # [ms**-1]

# (G_stretch, alpha_stretch, trans_p_mmHg, sigma0, E_SAC are included above in 
#  SMC flux Constants)

L_j = 0.025e-3                                   # [ms**-1] 

G_cat_j = 6.6e-7                                 # [uM mV**-1 ms**-1]
E_Ca_j = 50 #mV
m_3_cat_j = -0.18 
m_4_cat_j = 0.37 

G_tot_j = 6927e9                                 # [mpS] = [kg**-1*m**-2*ms**3*A**2] --> WRONG?
v_K_j = -80 #mV

c = -0.4
bb_j = -80.8 #mV
a_1_j = 53.3 #mV
a_2_j = 53.3 # mV 
m_3b_j = 1.32e-3 #mV**-1
m_4b_j = 0.3 #mV
m_3s_j = -0.28 
m_4s_j = 0.389 

G_R_j = 955e9                                    # [mpS] = [kg**-1*m**-2*ms**3*A**2] --> WRONG?
v_rest_j = -31.1 #mV
k_d_j = 0.1e-3                                   # [ms**-1]

P_Ca = 0.05e-3                                   # [ms**-1]
P_IP3 = 0.05e-3                                  # [ms**-1]
G_coup = 0.5e-3                                  # [ms**-1]

# Additional Equations Constants
alpha_act_i = 0.13 #uM**2
v_Ca3_i = -27 #mV
R_K_i = 12 #mV
z_1 = 4.5e-3 #mV uM**-1
z_2 = 112 #mV
z_3 = 4.2e-4 #uM**-1
z_4 = 12.6 # -
z_5 = -7.4e-2 #mV**-1

# NO pathway
K_mArg_j = 1.5 # [uM] 
K_mO2_j = 7.7 # [uM]  Chen2006
k_dno = 0.01e-3                                  # [ms**-1] 
K_m_mlcp = 5.5 # [uM] 
V_NOj_max = 1.22e-3                              # [ms**-1]  
LArg_j = 100 # [uM] 
k_O2 = 9.6e-9                                    # [uM**-2 ms**-1] 
W_0 = 1.4 #  shear gating constant (Comerford2008)
delta_wss = 2.86 # [Pa]  the membrane shear modulus (Comerford2008)
k_1 = 100e-3                                     # [ms**{-1}] 
k1 = 2                                           # [uM**-1 ms**-1] 
k2 = 0.1e-3                                      # [ms**-1] 
k3 = 3e-3                                        # [uM**-1 ms**-1] 
V_max_sGC = 0.8520e-3                            # [uM/ms] 
k_pde = 0.0195e-3                                # [ms**-1] 
C_4 = 0.011e-3                                   # [ms**{-1} microM**{-2}] 
K_m_pde = 2 # [uM] 
gam_eNOS = 0.1 # [-] 
mu2_j = 0.0167e-3                                # [ms**-1] 
K_dis = 9e-5                                     # [uM ms**-1] 
K_eNOS = 4.5e-1 # [uM]     
g_max = 0.06e-3                                  # [uM ms**-1]     
alp = 2 # [-]  zero shear open channel constant (Comerford2008 in Wiesner1997: alp = 3
delta_p_L = 9.1e-2 # Pa/um ME
x_ij = 3.75 # [um]  (Kavdia2002)

''' Wall mechanics '''
wallMech = 8.7                                   # scaling factor for the wall mechanics, original = 2.7 

# Contraction Equation Constants
K_3 = 0.4e-3                                     # [ms**-1]
K_4 = 0.1e-3                                     # [ms**-1]
K_7 = 0.1e-3                                     # [ms**-1]
gamma_cross = 17e-3                              # [uM**-3 ms**-1]
n_cross = 3                                      # [-]

# Mechanical Equation Constants
eta_R = 1e7                                      # [Pa ms]
R_init = 20                                      # [um] 
trans_p = 4000                                   # [Pa]
E_passive = 66e3                                 # [Pa]
E_active = 233e3                                 # [Pa]
alpha = 0.6                                      # [-]
k_mlcp_b = 0.0086e-3                             # [ms**-1]
k_mlcp_c = 0.0327e-3                             # [ms**-1]

''' 20-HETE parameters ''' #*** we should probably adjust these so 20-HETE doesn't have such a huge effect

AA_max = 29  # uM
AA_m = 0.161  # uM
Ca0 = 0.1432  # uM
D_AA = .5*.033    # um**2/ms, model estimate
tau_AA = x_ki ** 2 /  (2 * D_AA) # Time constant for AA diffusion between astrocyte and SMC

V_a = 0.212e-2  # ms**-1  old value 0.212e-3
K_a = 228.2  # uM
V_f = 0.0319e-2  # ms**-1 old value 0.0319e-3
K_f = 23.5  # uM
lambda_h = 2.0e-3   # ms**-1  old value 0.139e-3

H0 = 0.068   # 0.126 was wrong baseline value before, must have changed when other parameters changed

NO_rest = 0.02047
R_NO = 0.02

''' Oxygen extraction fraction parameters '''

f_in0 = 1.521           # Baseline f_in for normalising, note that CBF/CBF_init is not exactly 1 so we need to do this to get f_in = 1 at baseline #1.681
O2_0 = 0.01             # previously 0.02 (from Chang2013) [mM] baseline tissue O2 concentration
g_OEF = 0.2             # Ratio of tissue O2 to plasma O2 (0.01 / 0.053)
E_0 = 0.4               # [-] baseline OEF
J_0 = 0.053             # previously 0.032 [mM/s] steady state change in pump & background flux due to CBF

''' GABA parameters '''
GABAswitch = 1          # 1: GABA (opens chlorine channels on AC and SMC) and no neuronal activity

G_Tmin = 1              # [-] Min value of G_Tact (GABA-T activity), occurs when NO concentration is normal
G_Tmax = 2              # [-] Max value of G_Tact (GABA-T activity), occurs when NO = 0
GT_midpoint = 0.007     # [uM] midpoint of NO dependent G_Tact sigmoidal, M.E.
GT_slope = 0.003        # [uM] slope of NO dependent G_Tact sigmoidal, M.E.

g_slope = 0.15          # [-] Slope of GABA sigmoidal, M.E.
g_midpoint = 0.8        # [-] Midpoint of GABA sigmoidal, M.E.  
E_GABA = -75            # [mV] Reversal potential of GABA activated Cl channel
G_GABA = .24 * G_Cl_i    # Maximum conductance of GABA activated Cl channel, based on SMC Cl leak channel conductance, M.E.

beta_Glu = 4.2e-3       # clearance of Glu, M.E. 
beta_GABA = 4.2e-3      # clearance of GABA, M.E. (Maite)
GABAbase = 0            # baseline GABA concentration (Maite)

''' NPY parameters '''
NPYswitch = 1           # 1: NPY release (opens VOCCs and causes constriction)  

npy_increase = 0.06     # [-] proportion increase of VOCC conductance from baseline value due to NPY, e.g. 0.06 means 6% increase, M.E. 
npy_slope = 0.15        # [-] Slope of NPY sigmoidal, M.E. 
npy_midpoint = 0.8      # [-] Midpoint of NPY sigmoidal, M.E.

beta_NPY = 4.2e-3       # clearance of NPY, M.E. (Maite)
NPYbase = 0             # baseline NPY concentration (Maite)