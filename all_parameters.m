%% Parameters of the NVU model

function p = all_parameters()
p.num_variables = 51; % Number of state variables

%% Input [ms]
p.startpulse = 300e3;        % Start time of input stimulation
p.lengthpulse = 2e3;       % Length of stimulation
p.Tend = 330e3;             % End of simulation
p.dt = .01e3;               % Time step (overwritten to be smaller if p.Whiskerpad = 2)

%% Commonly changed parameters - move any parameters needed here (don't duplicate)

p.gamma_switch = 0;     % Turn on/off use of Gamma function. Probably shouldn't be used if lengthpulse is less than 10 or so
p.GluSwitch = 1;        % Turn on/off glutamate input
p.O2switch = 1;         % O2switch = 0 ATP is plentiful, O2switch = 1 ATP is limited (oxygen-limited regime)
p.NOswitch_NE = 0;      % Turn on/off neuronal NO production
p.NOswitch_EC_WSS = 0;  % Turn on/off WSS induced eNOS production
p.NOswitch_EC_CA = 0;   % Turn on/off Ca2+ induced eNOS production
p.wallMech = 2.7;       % scaling factor for the wall mechanics
p.LCpathway = 0;        % Locus coeruleus pathway used for Whiskerpad = 1,2,4
p.E_switch = 1;         % E_switch = 0: use Buxton1997 model for OEF (no dependency on tissue O2 conc.), E_switch = 1: OEF dependent on tissue O2 (currently not working)***

%% Whiskerpad: different cases
p.Whiskerpad = 0;       % 0: Model with square input pulse
                        % 1: Model with whiskerpad Zheng input interpolation
                        % 2: Model with whisker model with thalamic input (Pinto et al) with single input block (no trailing 2 sec stimulus as in Zheng) scaled by Zheng input
                        % 3: Model with 2 square input pulses for Zheng paradigm
                        % 4: Model with whiskerpad Zheng input interpolation and fitted parameters to experiment
                        % 5: model with inhibitory neurons only
                        
p.double_pulse = 0;     % 1: use the second pulse in the Zheng input data, 0: only one pulse. Used for p.Whiskerpad = 1 or 4

if p.Whiskerpad == 5 
    p.Pinput = 0.0; % inhibitory input only, no excitation
    p.Qinput = 1.0;
else
    p.Pinput = 1.0; % excitatory input only, no inhibition
    p.Qinput = 0.0;
end

if p.Whiskerpad == 1 || p.Whiskerpad == 2 || p.Whiskerpad == 3 || p.Whiskerpad == 4
    p.ISI = 7;                    % INDEX for time period between stimulations [0.6,1,2,3,4,6,8]
    ISI_array = [0.6,1,2,3,4,6,8]*1e3;    
    p.actual_ISI = ISI_array(p.ISI);            % time period between stimulations
    
    % p.stim = INDEX for length of initial stimulation [2,8,16] 
    if p.lengthpulse == 2e3
        p.stim = 1;
    elseif p.lengthpulse == 8e3
        p.stim = 2;
    elseif p.lengthpulse == 16e3
        p.stim = 3;
    end

end

if p.Whiskerpad == 1 || p.Whiskerpad == 3 || p.Whiskerpad == 4
    p.secondpulse = p.startpulse + p.lengthpulse + p.actual_ISI;  
    p.secondlength = 2e3;                       % length of second pulse
else
    p.actual_ISI = 0;       % No second pulse
    p.secondlength = 0;
end

if p.Whiskerpad == 2
    p.dt = 0.001e3; % Thalamic input with fast triangles needs a smaller timestep
end
  
% Alpha, beta and rest values adjusted for the different inputs
if p.Whiskerpad == 0 || p.Whiskerpad == 3 || p.Whiskerpad == 5
    p.alpha_K_e = 3.2;
    p.beta_K_e = 4.2e-3;
    p.alpha_Na_sa = 4.23;
    p.beta_Na_sa = 0.39e-3;
    p.alpha_Na_d = -2.12;
    p.beta_Na_d = 0.75e-3;
    p.K_eBase = 3.5;
    p.Na_saBase = 9.37;
    p.Na_dBase = 9.42;
    p.tau_e = 3;        
    p.tau_i = 3;           
elseif p.Whiskerpad == 1
    p.alpha_K_e = 3.5;
    p.beta_K_e = 10e-3;
    p.alpha_Na_sa = 4.23;
    p.beta_Na_sa = 0.65e-3;
    p.alpha_Na_d = -2.12;
    p.beta_Na_d = 0.7e-3;
    p.K_eBase = 3.567;
    p.Na_saBase = 9.237;
    p.Na_dBase = 9.318;
    p.tau_e = 3;        
    p.tau_i = 3;    
elseif p.Whiskerpad == 2 % Thalamic input
    p.alpha_K_e = 3*3.2;
    p.beta_K_e = 4.2e-3;
    p.alpha_Na_sa = 4.23;
    p.beta_Na_sa = 0.39e-3;
    p.alpha_Na_d = -2.12;
    p.beta_Na_d = 0.75e-3;
    p.K_eBase = 3.5;
    p.Na_saBase = 9.37;
    p.Na_dBase = 9.42;
    p.tau_e = 5;        
    p.tau_i = 15;    
elseif p.Whiskerpad == 4
    p.alpha_K_e = 3.3;    % 3.3 with wallMech = 2.0 for CBF fit to Zheng, 3.0 with wallMech = 2.0 for fit to Berwick
    p.beta_K_e = 10e-3;
    p.alpha_Na_sa = 4.23;
    p.beta_Na_sa = 0.65e-3;
    p.alpha_Na_d = -2.12;
    p.beta_Na_d = 0.7e-3;
    p.K_eBase = 3.567;
    p.Na_saBase = 9.237;
    p.Na_dBase = 9.318;
    p.tau_e = 3;        
    p.tau_i = 3;  
end

%% Neuron


if p.Whiskerpad == 2 % Thalamic input
    p.EI_relative = 0.5167; % [-]    Maximum value for abs(E-I) after external inputs
    p.EImin = 0.1396;       % [-]    Minimum value for abs(E-I) after external inputs
else
    p.EI_relative = 0.2894; % [-]    Maximum value for abs(E-I) after external inputs
    p.EImin = 0;       % [-]    Minimum value for abs(E-I) after external inputs
end

p.Tinput = 0.2;

% Gamma parameters
p.gamma_width = 6e3;   % [ms]
p.gamma_delay = -2e3;   % [ms]?

% Wilson and Cowan model parameters
% K_e, Na_sa, Na_d parameters
% p.K_eBase = 3.5;        % [mM]
% p.Na_saBase = 9.37;     % [mM]
% p.Na_dBase = 9.42;      % [mM]

% p.alpha_K_e = 3.2;      % [-]
% p.beta_K_e = 4.2e-3;    % [ms^-1]
% p.alpha_Na_sa = 4.23;   % [-]                     
% p.beta_Na_sa = 0.39e-3; % [ms^-1]
% p.alpha_Na_d = -2.12;   % [-]                    
% p.beta_Na_d = 0.75e-3;  % [ms^-1]

% E and I parameters
p.c1 = 12;              % [-]
p.c2 = 10;              % [-]
p.c3 = 13;              % [-]
p.c4 = 11;              % [-]
p.r_e = 1.0;            % [ms]?
p.r_i = 4.0;            % [ms]?
p.k_e_input = 10;
p.k_i_input = 10;
p.a_e = 1.2;            % [-]
p.theta_e = 2.8;        % [-]
p.a_i = 1.0;            % [-]
p.theta_i = 4.0;        % [-]

% Thalamic input parameters
p.wee = 42;
p.wie = 24.6;
p.wei = 42;
p.wii = 18;
p.wte = 53.43;
p.wti = 68.4;

% p.Pinput = 1.0;         % [-]
% p.Qinput = 0.0;         % [-]

% Glutamate parameters
p.Glu_max = 1846;       % [uM] (one vesicle, Santucci2008)
p.Glu_slope = 0.1;      % [mM] Slope of sigmoidal (M.E.)
p.Ke_switch = 5.5;      % [mM] Threshold past which glutamate vesicle is released (M.E.)                               

% Elshin model constants
p.k_syn = 11.5;         % [-] Coupling scaling factor, values can range between 1.06 to 14.95 according to estimation from experimental data of Ventura and Harris
p.Farad = 96.485;       % [C / mmol]
p.ph = 26.6995;         % [RT/F]

% BOLD constants
p.tau_MTT = 3e3;        % [ms]
p.tau_TAT = 20e3;       % [ms]
p.d = 0.4;              % [-]
p.a_1 = 3.4;            % [-]   
p.a_2 = 1;              % [-]
p.V_0 = 0.03;           % [-]

% Oxygen model
p.alpha_O2 =  0.05;     % [-] Percentage of ATP production independent of O2 0.05;
p.K_init_e = 2.9;       % [mM]
p.Na_init_sa = 10;      % [mM]
p.Na_init_d = 10;       % [mM]

p.CBF_init = 3.2e-2;    % [-]3.2e-2
p.gamma_O2 = 0.1;       % [-] 0.1 

% NO pathway 
p.m_c = 4;              % [-] 
p.K_mA = 650;           % [uM] - fit to Santucci2008
p.K_mB = 2800;          % [uM] - fit to Santucci2008
p.v_n = -40;            % [mV] ; the neuronal membrane potential , assumed to be approx constant in this model
p.G_M = 0.46;           % original value 46e9 !! [pS] = [kg^-1*m^-2*ms^3*A^2] --> converted to ms
p.P_Ca_P_M = 3.6;       % [-] ; the relative conductance of the NMDA channel to Ca2+ compared to monovalent ions
p.Ca_ex = 2e3;          % [uM] ; the external calcium concentration (in Comerford+David2008: 1.5 mM!)
p.M = 1.3e5;            % [uM] ; the concentration of monovalent ions in the neuron
p.n_NR2A = 0.63;        % [-] ; average number of NR2A NMDA receptors per synapse (Santucci2008)
p.n_NR2B = 11;          % [-] ; average number of NR2B NMDA receptors per synapse (Santucci2008)
p.V_max_NO_n = 4.22e-3; % [ms^-1] ; maximum catalytic rate of NO production (Chen2006) - obtained from fig 6 & equ 17 & 18
p.O2_n = 200;           % [uM] ; tissue O2 concentration in the neuron (M.E.)
p.K_mO2_n = 243;        % [uM] ; Chen2006
p.LArg_n = 100;         % [uM] ; 
p.K_mArg_n = 1.5;       % [uM] ;     
p.k_O2_n = 9.6e-9;      % [uM^-2 ms^-1] ; % (Kavdia2002)
p.x_nk = 25;            % [um] ;  (M.E.)
p.D_cNO = 3300e-3;      % [um^2 ms^-1] ; Diffusion coefficient NO (Malinski1993)   
p.V_spine = 8e-5;       % [pL] ; volume of the neuronal dendritic spine Santucci2008
p.k_ex = 1600e-3;       % [ms^-1] ; decay rate constant of internal calcium concentration Santucci2008
p.Ca_rest = 0.1;        % [uM] ; resting calcium concentration (in Comerford+David2008: 2.830 mM; in Santucci2008P: 0.1 \muM)
p.lambda_buf = 20;      % [-] ; buffer capacity Santucci2008
p.V_maxNOS = 0.697e-3;     % [uM ms^-1] ; M.E. original values 25e-6
p.K_actNOS = 0.8;   % [uM] ; original values 9.27e-5
p.mu2_n = 0.2e-2;    % [ms^-1] ; original values 0.0167e-3 

%% Astrocyte
    
% Conductances - now in uM mV^-1 ms^-1 so consistent with other
% fluxes, original conductances in mho m^-2 are commented
% Converted by dividing by R_k=6e-8 and F=9.65e4
p.G_BK_k = 10.25e-3;        % g_BK_k = 225 * 1e-12 / 3.7e-9 = 6.08e-2;
p.G_K_k = 6907.77e-3;      % 40
p.G_Na_k = 226.94e-3;      % 1.314
p.G_NBC_k = 130.74e-3;     % 7.57e-1
p.G_KCC1_k = 1.728e-3;     % 1e-2
p.G_NKCC1_k = 9.568e-3;    % 5.54e-2
p.G_Cl_k = 151.93e-3;      % 8.797e-1

p.J_NaK_max = 2.3667e1;                             % [uM ms^-1]      1.42e-3

% Switches to turn on and off some things
p.rhoSwitch = 1; 
p.trpv_switch = 1; 

p.rho_min = 0.1;     
p.rho_max = 0.7;  

p.v_4 = 8; %mV
p.v_5 = 15; %mV
p.v_6 = -55; %mV
p.Ca_3 = 0.4; % uM 
p.Ca_4 = 0.35; % uM 

% Scaling Constants
p.R_s = 2.79e-8; % m
p.R_k = 6e-8; % m
p.R_tot = 8.79e-8; % m
p.VR_sa = p.R_s/p.R_k;

% Calcium in the Astrocyte Equations Constants
p.delta = 1.235e-2; 
p.VR_ER_cyt = 0.185;
p.k_on = 2e-3;                                      % [uM ms^-1]
p.K_inh = 0.1; %uM
p.r_h = 4.8e-3;                                     % [uM/ms]
p.k_deg = 1.25e-3;                                  % [ms^-1]
p.V_eet = 72e-3;                                    % [ms^-1]
p.k_eet = 7.2e-3;                                   % [ms^-1]
p.Ca_k_min = 0.1; % uM
p.eet_shift = 2; %mV/uM
p.K_I = 0.03; % uM

%TRPV4
p.Capmin_k = 2000; %uM
p.C_astr_k = 40;%pF
p.gam_cae_k = 200; %uM
p.gam_cai_k = 0.01; %uM
p.epshalf_k = 0.1; % 
p.kappa_k = 0.1;
p.v1_TRPV_k = 120; %mV
p.v2_TRPV_k = 13; %mV
p.t_TRPV_k = 0.9e3;                                 % [ms] 
p.Ca_decay_k = 0.5e-3;                              % [ms^-1]
p.G_TRPV_k = 3.15e-7;                               % [ms^-1]
p.r_buff = 0.05; % Rate at which Ca2+ from the TRPV4 channel at the endfoot is buffered compared to rest of channels on the astrocyte body [-]

% Perivascular space
p.VR_pa = 0.001;% [-]
p.VR_ps = 0.001;% [-]
p.R_decay = 0.15e-3;                                % [ms^-1]
p.K_p_min = 3e3;% uM

% Fluxes Constants
p.R_g = 8.315; %J mol^-1 K^-1; Gas constant
p.T = 300; % K; Temperature

p.K_Na_k = 10000; % uM
p.K_K_s = 1500; % uM

p.J_max = 2880e-3;                                  % [uM ms^-1]
p.K_act = 0.17; %uM
p.P_L = 0.0804e-3;                                  % [uM ms^-1]
p.V_max = 20e-3;                                    % [uM ms^-1]
p.k_pump = 0.24; %uM

% Additional Equations; Astrocyte Constants
p.z_K = 1;% [-]
p.z_Na = 1;% [-]
p.z_Cl = -1;% [-]
p.z_NBC = -1;% [-]
p.z_Ca = 2;% [-]
p.BK_end = 40;% [-]
p.K_ex = 0.26; %uM
p.B_ex = 11.35; %uM
p.K_G = 8.82; 
p.psi_w = 2.664e-3;                                 % [ms^-1]

% NO Pathway
p.x_ki = 25;            % [um] ;  (M.E.)
p.k_O2_k = 9.6e-9;      % [uM^-2 ms^-1] ;  (Kavdia2002)
p.O2_k = 200;           % [uM] ;  (M.E.)

%% SMC/EC

% Smooth Muscle Cell ODE Constants
p.gamma_i = 1970; %mV uM^-1
p.lambda_i = 45e-3;                                 % [ms^-1]

% Endothelial Cell ODE Constants
p.C_m_j = 25.8;                                     % pF veranderen naar iets met ms??
p.J_PLC = 0.11e-3;%0.11 for steady state, 0.3 for oscillations % [uM ms^-1] 
p.J_0_j = 0.029e-3;                                 % [uM ms^-1] constant Ca influx (EC)

% Smooth Muscle Cell Flux Constants
p.F_i = 0.23e-3;                                    % [uM ms^-1]      % IP3/RYR channel strength
p.K_r_i = 1; %uM

p.B_i = 2.025e-3;                                   % [uM ms^-1]     % SERCA pump strength     
p.c_b_i = 1.0; %uM

p.C_i = 55e-3;                                      % [uM ms^-1]
p.s_c_i = 2.0; %uM
p.c_c_i = 0.9; %uM

p.D_i = 0.24e-3;                                    % [ms^-1]
p.v_d = -100; %mV
p.R_d_i = 250; %mV

p.L_i = 0.025e-3;                                   % [ms^-1]

p.G_Ca_i = 1.29e-6;                                 % [uM mV^-1 ms^-1]
p.v_Ca1_i = 100; %mV
p.v_Ca2_i = -24; %mV
p.R_Ca_i = 8.5; %mV

p.G_NaCa_i = 3.16e-6;                               % [uM mV^-1 ms^-1]
p.c_NaCa_i = 0.5; %uM
p.v_NaCa_i = -30; %mV

p.G_stretch = 6.1e-6;                               % [uM mV^-1 ms^-1]   (Also EC parameter)
p.alpha_stretch = 7.4e-3; % mmHg^-1     (Also EC parameter)
p.trans_p_mmHg = 30; % mmHg             (Also EC parameter) transmural pressure. 30 mmHg = 4000 Pa
p.sigma_0 = 500; % mmHg                 (Also EC parameter)
p.E_SAC = -18; % mV                     (Also EC parameter)

p.F_NaK_i = 4.32e-5;                                % [uM ms^-1]

p.G_Cl_i = 1.34e-6;                                 % [uM mV^-1 ms^-1]
p.v_Cl_i = -25; %mV

p.G_K_i = 4.46e-6;                                  % [uM mV^-1 ms^-1]
p.v_K_i = -94; %mV

p.F_KIR_i = 1.285e-9;%0.381;                       % [uM mV^-1 ms^-1] F_KIR_i = 1.285e-6
p.k_d_i = 0.1e-3;                                   % [ms^-1]

% Endothelial Cell Flux Constants
p.F_j = 0.23e-3;                                    % [uM ms^-1]
p.K_r_j = 1; %uM
p.B_j = 0.5e-3;                                     % [uM ms^-1]
p.c_b_j = 1; %uM
p.C_j = 5e-3;                                       % [uM ms^-1]
p.s_c_j = 2; %uM
p.c_c_j = 0.9; %uM
p.D_j = 0.24e-3;                                    % [ms^-1]

% (G_stretch, alpha_stretch, trans_p_mmHg, sigma0, E_SAC are included above in 
%  SMC flux Constants)

p.L_j = 0.025e-3;                                   % [ms^-1] 

p.G_cat_j = 6.6e-7;                                 % [uM mV^-1 ms^-1]
p.E_Ca_j = 50; %mV
p.m_3_cat_j = -0.18; 
p.m_4_cat_j = 0.37; 

p.G_tot_j = 6927e9;                                 % [mpS] = [kg^-1*m^-2*ms^3*A^2] --> WRONG?
p.v_K_j = -80; %mV

p.c = -0.4;
p.bb_j = -80.8; %mV
p.a_1_j = 53.3; %mV
p.a_2_j = 53.3; % mV 
p.m_3b_j = 1.32e-3; %mV^-1
p.m_4b_j = 0.3; %mV
p.m_3s_j = -0.28; 
p.m_4s_j = 0.389; 

p.G_R_j = 955e9;                                    % [mpS] = [kg^-1*m^-2*ms^3*A^2] --> WRONG?
p.v_rest_j = -31.1; %mV
p.k_d_j = 0.1e-3;                                   % [ms^-1]

p.P_Ca = 0.05e-3;                                   % [ms^-1]
p.P_IP3 = 0.05e-3;                                  % [ms^-1]
p.G_coup = 0.5e-3;                                  % [ms^-1]

% Additional Equations Constants
p.alpha_act_i = 0.13; %uM^2
p.v_Ca3_i = -27; %mV
p.R_K_i = 12; %mV
p.z_1 = 4.5e-3; %mV uM^-1
p.z_2 = 112; %mV
p.z_3 = 4.2e-4; %uM^-1
p.z_4 = 12.6; % -
p.z_5 = -7.4e-2; %mV^-1

% NO pathway
p.K_mArg_j = 1.5; % [uM] ;
p.K_mO2_j = 7.7; % [uM] ; Chen2006
p.k_dno = 0.01e-3;                                  % [ms^-1] ;
p.K_m_mlcp = 5.5; % [uM] ;
p.V_NOj_max = 1.22e-3;                              % [ms^-1] ; 
p.LArg_j = 100; % [uM] ;
p.k_O2 = 9.6e-9;                                    % [uM^-2 ms^-1] ;
p.W_0 = 1.4; %  shear gating constant (Comerford2008)
p.delta_wss = 2.86; % [Pa] ; the membrane shear modulus (Comerford2008)
p.k_1 = 100e-3;                                     % [ms^{-1}] ;
p.k1 = 2;                                           % [uM^-1 ms^-1] ;
p.k2 = 0.1e-3;                                      % [ms^-1] ;
p.k3 = 3e-3;                                        % [uM^-1 ms^-1] ;
p.V_max_sGC = 0.8520e-3;                            % [uM/ms] ;
p.k_pde = 0.0195e-3;                                % [ms^-1] ;
p.C_4 = 0.011e-3;                                   % [ms^{-1} microM^{-2}] ;
p.K_m_pde = 2; % [uM] ;
p.gam_eNOS = 0.1; % [-] ;
p.mu2_j = 0.0167e-3;                                % [ms^-1] ;
p.K_dis = 9e-5;                                     % [uM ms^-1] ;
p.K_eNOS = 4.5e-1; % [uM] ;    
p.g_max = 0.06e-3;                                  % [uM ms^-1] ;    
p.alp = 2; % [-] ; zero shear open channel constant (Comerford2008; in Wiesner1997: alp = 3
p.delta_p_L = 9.1e-2; % Pa/um ME
p.x_ij = 3.75; % [um]  (Kavdia2002)



%% Wall mechanics

% Contraction Equation Constants
p.K_3 = 0.4e-3;                                     % [ms^-1]
p.K_4 = 0.1e-3;                                     % [ms^-1]
p.K_7 = 0.1e-3;                                     % [ms^-1]
p.gamma_cross = 17e-3;                              % [uM^-3 ms^-1]
p.n_cross = 3;                                      % [-]

% Mechanical Equation Constants
p.eta_R = 1e7;                                        % [Pa ms]
p.R_init = 19;%20;                                      % [um]
p.trans_p = 4000;                                   % [Pa]
p.E_passive = 66e3;                                 % [Pa]
p.E_active = 233e3;                                 % [Pa]
p.alpha = 0.6;                                      % [-]
p.k_mlcp_b = 0.0086e-3;                             % [ms^-1]
p.k_mlcp_c = 0.0327e-3;                             % [ms^-1]

%% 20-HETE parameters
p.AA_max = 29;  % uM
p.AA_m = 0.161;  % uM
p.Ca0 = 0.1432;  % uM
p.D_AA = 0.033;    % um^2/ms, model estimate
p.tau_AA = p.x_ki ^ 2 ./  (2 * p.D_AA); % Time constant for AA diffusion between astrocyte and SMC

p.V_a = 0.212e-2;  % ms^-1  old value 0.212e-3
p.K_a = 228.2;  % uM
p.V_f = 0.0319e-2;  % ms^-1 old value 0.0319e-3
p.K_f = 23.5;  % uM
p.lambda_h = 2.0e-3;  % ms^-1  old value 0.139e-3
p.Hshift = 30;
p.H0 = 0.126;
p.NO_rest = 0.02047;
p.R_NO = 0.02;

%% Oxygen extraction fraction parameters
% p.O2_b = 5.175;             % previously 4e-2 (from Chang2013); [mM] baseline capillary O2 concentration (not used anymore)
p.O2_0 = 0.01;              % previously 0.02 (from Chang2013); [mM] baseline tissue O2 concentration
p.g_OEF = 0.2;              % Ratio of tissue O2 to plasma O2 (0.01 / 0.053)
p.E_0 = 0.4;                % [-] baseline OEF
p.J_0 = 0.053;             % previously 0.032; [mM/s] steady state change in pump & background flux due to CBF
end
