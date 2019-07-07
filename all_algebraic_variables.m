%% Fluxes for the NVU system
% Make sure to use .*, ./, and .^ instead of *, /, ^ (multiplication, division and exponentiation)

%% Neuron

if p.Whiskerpad == 2 % Thalamic input
    
    T = T(t,p,input_data)';
    arg_e = p.wee .* E_t - p.wie .* I_t + p.wte .* T; 
    arg_i = p.wei .* E_t - p.wii .* I_t + p.wti .* T; 
    f_arg_e = fe(p,arg_e);
    f_arg_i = fi(p,arg_i);
    
else % Normal P,Q input
    
    k_e = se(p,p.k_e_input); %1 ./ (1 + exp(-(p.a_e .* (p.k_e_input - p.theta_e)))) - 1 ./ (1 + exp(p.a_e .* p.theta_e));
    k_i = si(p,p.k_e_input); %1 ./ (1 + exp(-(p.a_i .* (p.k_i_input - p.theta_i)))) - 1 ./ (1 + exp(p.a_i .* p.theta_i));

    P = P(t,p,input_data);
    Q = Q(t,p,input_data);
    if size(E_t) ~= size(P) % Minor check to see if they are the same dimensions (if not then invert P/Q), only an issue when trying to plot elements as P/Q swap dimensions when using input_data for some reason
        P = P';  % [-]           Excitatory cells response function input
        Q = Q';  % [-]           Inhibitory cells response function input
    end

    arg_e = p.c1 .* E_t - p.c2 .* I_t + P;  % [-]           Excitatory cells response function input
    arg_i = p.c3 .* E_t - p.c4 .* I_t + Q;  % [-]           Inhibitory cells response function input

    s_arg_e = se(p,arg_e); %1 ./ (1 + exp(-(p.a_e .* (arg_e - p.theta_e)))) - 1 ./ (1 + exp(p.a_e .* p.theta_e));
    s_arg_i = si(p,arg_i); %1 ./ (1 + exp(-(p.a_i .* (arg_i - p.theta_i)))) - 1 ./ (1 + exp(p.a_i .* p.theta_i));
end

O2_p            = p.O2_0 * (1 - p.O2switch) + O2 * p.O2switch;  % [mM]  Oxygen concentration dependent on ATP                                   
J_pump2         = 2 * (1 + p.O2_0 ./ (((1 - p.alpha_O2) * O2_p)  + p.alpha_O2 * p.O2_0)).^(-1); % [mM s^-1]       
J_pump2_0       = 0.0952;           % [mM s^-1]                      
J_pump2_O2_0    = 1;                % [mM s^-1]                         
P_02            = (J_pump2 - J_pump2_0 ) ./ ( J_pump2_O2_0 - J_pump2_0);   % [-]  Normalised pump rate of oxygen            

CBF             = p.CBF_init * (R.^4 / p.R_init^4);             % [-] Normalised cerebral blood flow
% CBF             = p.CBF_init * (R.^4 / 22.41^4); % Actual baseline value so CBF/CBF_init = 1      
                            
J_O2_background = p.J_0 * P_02 * (1 - p.gamma_O2);  % [mM s^-1]    Background oxygen consumption                                                                                 
J_pump1_sa      = (1 + (p.K_init_e ./ K_e)).^(-2) .* (1 +  (p.Na_init_sa ./ Na_sa)) .^ (-3); % [mM s^-1]  Oxygen consumption in soma/axon pump
J_pump1init_sa  = 0.0312;           % [mM s^-1]   Initial consumption in soma/axon pump    
J_pump1_d       = (1 + (p.K_init_e ./ K_e)).^(-2) .* (1 + (p.Na_init_d ./ Na_d)).^(-3);   % [mM s^-1]  Oxygen consumption in dendrite pump
J_pump1init_d   = 0.0312;           % [mM s^-1]   Initial consumption in dendrite pump             
J_O2_pump       = p.J_0 .* P_02 .* p.gamma_O2 .* ((J_pump1_sa + J_pump1_d) ./ (J_pump1init_sa + J_pump1init_d)); % [mM s^-1]    Total oxygen consumption of the NaK pumps  

% BOLD
f_out           = CBV.^(1/p.d) + p.tau_TAT * (1/(p.tau_MTT + p.tau_TAT) .* ( CBF/p.CBF_init  - CBV.^(1/p.d) ));  % [-]   Buxton balloon model outflow                                     
% CMRO2           = J_O2_background + J_O2_pump;  % [mM s^-1]   Metabolic rate of oxygen in the neurons                                 
% CMRO2_init      = p.CBF_init * P_02;% [mM s^-1]  Initial metabolic rate of oxygen

% NO pathway
Glu = p.GluSwitch * 0.5 * p.Glu_max * ( 1 + tanh( (K_e - p.Ke_switch) / p.Glu_slope) );
w_NR2A = Glu ./ (p.K_mA + Glu);     % [-] 
w_NR2B = Glu ./ (p.K_mB + Glu);     % [-] 
I_Ca = (-4 * p.v_n * p.G_M * p.P_Ca_P_M * (p.Ca_ex / p.M)) / (1 + exp(-0.08 * (p.v_n + 20))) * (exp(2 * p.v_n / p.ph)) / (1 - exp(2 * p.v_n / p.ph));   % [fA]                                    
I_Ca_tot = I_Ca .* (p.n_NR2A * w_NR2A + p.n_NR2B * w_NR2B);   % [fA]                                                 
CaM = Ca_n / p.m_c;                 % [uM]        
p_NO_n = nNOS_act_n * p.V_max_NO_n * p.O2_n /  (p.K_mO2_n + p.O2_n) * p.LArg_n / (p.K_mArg_n + p.LArg_n);
c_NO_n = p.k_O2_n * NO_n.^2 * p.O2_n;   
tau_nk = p.x_nk ^ 2 ./  (2 * p.D_cNO);                          % [ms]
d_NO_n = (NO_k - NO_n) ./ tau_nk;         

% Flux from ECS to SC
if p.gamma_switch == 0
    dKedt = -p.beta_K_e * (K_e - p.K_eBase) + (p.alpha_K_e .* p.beta_K_e) .* ((abs(E_t - I_t) - p.EImin)./ (p.EI_relative - p.EImin));          % K_e derivative as neuron activation in Ostby            
else
    dKedt = - p.beta_K_e * (K_e - p.K_eBase) + Gamma(t,p,E_t,I_t); 
end

J_K_NEtoSC = p.k_syn * dKedt * 1e3;                                   % Input flux: K+ flux into SC
J_Na_NEtoSC = -p.k_syn * dKedt * 1e3;                                 % Input flux: Na+ flux out of SC

%% Astrocyte

% Electroneutrality condition
Cl_s = Na_s + K_s - HCO3_s;

% Nernst potentials (in mV)
E_K_k = p.ph / p.z_K * log(K_s ./ K_k);
E_Na_k = p.ph / p.z_Na * log(Na_s ./ Na_k);
E_Cl_k = p.ph / p.z_Cl * log(Cl_s ./ Cl_k);
E_NBC_k = p.ph / p.z_NBC * log((Na_s .* HCO3_s.^2) ./ (Na_k .* HCO3_k.^2));
E_BK_k = p.ph / p.z_K * log(K_p ./ K_k);        
E_TRPV_k = p.ph / p.z_Ca * log(Ca_p ./ Ca_k);                   % Nernst potential TRPV

% Flux through the Sodium Potassium pump
J_NaK_k = p.J_NaK_max * Na_k.^1.5 ./ (Na_k.^1.5 + p.K_Na_k^1.5) .* K_s ./ (K_s + p.K_K_s);

% Fluxes
J_BK_k = p.G_BK_k * w_k .* (v_k - E_BK_k);                      % BK flux (uM/ms)
J_BK_p = J_BK_k ./ p.VR_pa;                                     % K+ influx into the PVS (uM/ms)
J_K_k = p.G_K_k * (v_k - E_K_k);
J_Na_k = p.G_Na_k * (v_k - E_Na_k);
J_NBC_k = p.G_NBC_k * (v_k - E_NBC_k);
J_KCC1_k = p.G_KCC1_k * p.ph .* log((K_s .* Cl_s) ./ (K_k .* Cl_k));
J_NKCC1_k = p.G_NKCC1_k * p.ph .* log((Na_s .* K_s .* Cl_s.^2) ./ (Na_k .* K_k .* Cl_k.^2));
J_Cl_k = p.G_Cl_k * (v_k - E_Cl_k);

% Calcium Equations
% Flux
J_IP3 = p.J_max * ( I_k ./ (I_k + p.K_I) .*  Ca_k ./ (Ca_k + p.K_act) .* h_k).^3 .* (1 - Ca_k ./ s_k);
J_ER_leak = p.P_L * (1 - Ca_k ./ s_k);
J_pump = p.V_max * Ca_k.^2 ./ (Ca_k.^2 + p.k_pump^2);
J_TRPV_k = p.G_TRPV_k * m_k .* (v_k - E_TRPV_k);

rho = p.rho_min + (p.rho_max - p.rho_min)/p.Glu_max * Glu;

% Other equations
B_cyt = 1 ./ (1 + p.BK_end + p.K_ex * p.B_ex ./ (p.K_ex + Ca_k).^2);
G = (rho + p.delta) ./ (p.K_G + rho + p.delta);
v_3 = p.v_6 - p.v_5 / 2 * tanh((Ca_k - p.Ca_3) / p.Ca_4);

% Parent Calcium equations 
w_inf = 0.5 * (1 + tanh((v_k + p.eet_shift * eet_k - v_3) / p.v_4));
phi_w = p.psi_w * cosh((v_k - v_3) / (2 * p.v_4));

% TRPV Channel open probability equations
H_Ca_k = Ca_k ./ p.gam_cai_k + Ca_p ./ p.gam_cae_k;
eta = (R - p.R_init) ./ (p.R_init);
minf_k = (1 ./ (1 + exp(-(eta - p.epshalf_k) ./ p.kappa_k))) .* ((1 ./ (1 + H_Ca_k)) .* (H_Ca_k + tanh((v_k - p.v1_TRPV_k) ./ p.v2_TRPV_k))); 

% NO pathway
tau_ki = p.x_ki ^ 2 ./  (2 * p.D_cNO);                          % [ms]
p_NO_k = 0;
c_NO_k = p.k_O2_k * NO_k.^2 * p.O2_k;                           % [uM/ms]
d_NO_k = (NO_n - NO_k) ./ tau_nk + (NO_i - NO_k) ./ tau_ki;     % [uM/ms]

%% SMC/EC

% SMC fluxes
J_IP3_i = p.F_i * I_i.^2 ./ (p.K_r_i^2 + I_i.^2);           % IP3 channel
J_SR_uptake_i = p.B_i * Ca_i.^2 ./ (p.c_b_i^2 + Ca_i.^2);   % SERCA pump
J_CICR_i = p.C_i * s_i.^2 ./ (p.s_c_i^2 + s_i.^2) .* Ca_i.^4 ./ (p.c_c_i^4 + Ca_i.^4);
J_extrusion_i = p.D_i * Ca_i .* (1 + (v_i - p.v_d) / p.R_d_i);
J_SR_leak_i = p.L_i * s_i;
J_VOCC_i = p.G_Ca_i .* (v_i - p.v_Ca1_i) ./ (1 + exp(-(v_i - p.v_Ca2_i) ./ p.R_Ca_i)); %***
J_NaCa_i = p.G_NaCa_i * Ca_i ./ (Ca_i + p.c_NaCa_i) .* (v_i - p.v_NaCa_i);

h = 0.1 * R; 
J_stretch_i = p.G_stretch ./ (1 + exp(-p.alpha_stretch*(p.trans_p_mmHg*R./h - p.sigma_0))) .* (v_i - p.E_SAC);

J_Cl_i = p.G_Cl_i * (v_i - p.v_Cl_i);
J_NaK_i = p.F_NaK_i;
J_K_i   = p.G_K_i * w_i .* (v_i - p.v_K_i);
J_degrad_i = p.k_d_i * I_i;

v_KIR_i = p.z_1 * K_p - p.z_2;
G_KIR_i = p.F_KIR_i * exp(p.z_5 * v_i + p.z_3 * K_p);           % Changed by including exp(-z_4) parameter into F_KIR_i parameter, no large constant in exponential equation
J_KIR_i = G_KIR_i .* (v_i - v_KIR_i);

% EC fluxes
J_IP3_j = p.F_j * I_j.^2 ./ (p.K_r_j^2 + I_j.^2);
J_ER_uptake_j = p.B_j * Ca_j.^2 ./ (p.c_b_j^2 + Ca_j.^2);  
J_CICR_j = p.C_j * s_j.^2 ./ (p.s_c_j^2 + s_j.^2) .* Ca_j.^4 ./ (p.c_c_j^4 + Ca_j.^4);
J_extrusion_j = p.D_j * Ca_j;
J_stretch_j = p.G_stretch ./ (1 + exp(-p.alpha_stretch*(p.trans_p_mmHg*R./h - p.sigma_0))) .* (v_j - p.E_SAC);
J_ER_leak_j = p.L_j * s_j;
J_cation_j = p.G_cat_j * (p.E_Ca_j - v_j) * 0.5 .* (1 + tanh((log10(Ca_j) - p.m_3_cat_j) / p.m_4_cat_j));
J_BK_Ca_j = 0.2 * (1 + tanh( ((log10(Ca_j) - p.c) .* (v_j - p.bb_j) - p.a_1_j) ./ (p.m_3b_j * (v_j + p.a_2_j*(log10(Ca_j) - p.c) - p.bb_j).^2 + p.m_4b_j)));
J_SK_Ca_j = 0.3 * (1 + tanh((log10(Ca_j) - p.m_3s_j) / p.m_4s_j));
J_K_j = p.G_tot_j * (v_j - p.v_K_j) .* (J_BK_Ca_j + J_SK_Ca_j);
J_R_j = p.G_R_j * (v_j - p.v_rest_j);
J_degrad_j = p.k_d_j * I_j;

% Coupling
V_coup_i = -p.G_coup * (v_i - v_j);
J_IP3_coup_i = -p.P_IP3 * (I_i - I_j);
J_Ca_coup_i = -p.P_Ca * (Ca_i - Ca_j);

c_w_i = 1/2 * (1 + tanh((cGMP_i - 10.75)/0.668) );
K_act_i = (Ca_i + c_w_i).^2 ./ ((Ca_i + c_w_i).^2 + p.alpha_act_i * exp(-(v_i - p.v_Ca3_i - p.Hshift*(H_i-p.H0)) / p.R_K_i)); % Is shifted by 20-HETE
tau_wss = R/2 * p.delta_p_L; % from Dormanns 2016

% NO pathway 
tau_ki = p.x_ki ^ 2 ./  (2 * p.D_cNO);                          % [ms]
tau_ij = p.x_ij ^ 2 ./  (2 * p.D_cNO);                          % [ms]
p_NO_i = 0;
c_NO_i = p.k_dno * NO_i;
d_NO_i = (NO_k - NO_i) ./ tau_ki + (NO_j - NO_i) ./ tau_ij;     % [uM]
k4 = p.C_4 * cGMP_i.^2;
E_5c = 1 - E_b - E_6c;
V_max_pde = p.k_pde * cGMP_i;
R_cGMP2 = cGMP_i.^2 ./ (cGMP_i.^2 + p.K_m_mlcp^2);

O2_j = O2*1e3;  % Oxygen in EC taken as O2 from lumen (diffusion very fast so plausible!) instead of constant, in uM
p_NO_j = ( p.V_NOj_max .* eNOS_act_j .* O2_j ./ (p.K_mO2_j + O2_j) .* p.LArg_j ./ (p.K_mArg_j + p.LArg_j) );
c_NO_j = p.k_O2 .* NO_j.^2 .* O2_j;% consumption by oxygen 

J_lumen = - NO_j * 4 * p.D_cNO ./ (25.^2); 
d_NO_j = (NO_i - NO_j) ./ tau_ij + J_lumen; 

W_wss = p.W_0 * (tau_wss + sqrt(16 * p.delta_wss^2 + tau_wss.^2) - 4 * p.delta_wss).^2 ./ (tau_wss + sqrt(16 * p.delta_wss.^2 + tau_wss.^2)); 
F_wss = 1 ./ (1 + p.alp .* exp(- W_wss)) - 1 ./ (1 + p.alp); 
Act_eNOS_Ca = p.K_dis .* Ca_j ./ (p.K_eNOS + Ca_j); 
Act_eNOS_wss = p.g_max .* F_wss;  

% 20-HETE
f_NO = 1 ./ (1 + exp((NO_i - p.NO_rest) ./ p.R_NO));

%% Wall mechanics

K_1 = p.gamma_cross * Ca_i.^p.n_cross;                          
K_6 = K_1;                                                      
K_2 = 58.1395 * p.k_mlcp_b + 58.1395 * p.k_mlcp_c * R_cGMP2;    
K_5 = K_2;                                                      

M = 1 - AM - AMp - Mp;        

% Mechanical Equations
F_r = AMp + AM;
E = p.E_passive + F_r * (p.E_active - p.E_passive);
R_0 = p.R_init + F_r * (p.alpha - 1) * p.R_init;

%% Oxygen extraction fraction & CMRO2 equations
 
f_in_dim = CBF./(p.CBF_init);                                                      % inflow of blood (normalised CBF)

f_in = f_in_dim./p.f_in0;

if p.E_switch == 0  
    OEF = 1 - (1 - p.E_0).^(1./f_in);                                       % Use simplified model of OEF based on Buxton 1997 (no dependency on tissue O2)
else
    E1 = 1-(1-p.E_0).^(1-p.g_OEF);
    OEF = 1 - (1 - E1).^(1./(f_in.*(1 - p.g_OEF)));                         % Oxygen extraction fraction dependent on g = [tissue O2]/[plasma O2]
end

CMRO2 = f_in .* OEF ./ p.E_0;                                               % Cerebral metabolic rate of oxygen consumption
J_O2_vascular   = CBF .* OEF ./ p.E_0;                                      % [mM s^-1] Vascular supply of oxygen, previously CBF .* ((p.O2_b - O2) ./ (p.O2_b - p.O2_0));   
