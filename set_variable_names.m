%% Set variable names for plotting in run script

% Neuron
E_t =   u(:, idx.E_t);              % [-]                       Proportion of excitatory cells firing
I_t =   u(:, idx.I_t);              % [-]                       Proportion of inhibitory cells firing
K_e =   u(:, idx.K_e);              % [mM]                      K+ concentration of ECS
Na_sa = u(:, idx.Na_sa);            % [mM]                      Na+ concentration of soma/axon
Na_d =  u(:, idx.Na_d);             % [mM]                      Na+ concentration of dendrite 
O2 =    u(:, idx.O2);               % [mM]                      Tissue oxygen concentration
CBV =   u(:, idx.CBV);              % [-]                       BOLD CBV
HbR =   u(:, idx.HbR);              % [-]                       BOLD deoxyhemoglobin
Ca_n =  u(:, idx.Ca_n);             % [uM]                      Calcium concentration     
nNOS_act_n = u(:, idx.nNOS_act_n);  % [uM]                      nNOS active concentration
NO_n =  u(:, idx.NO_n);             % [uM]                      NO concentration

% Astrocyte
v_k = u(:, idx.v_k);
K_p = u(:, idx.K_p);
Ca_p = u(:, idx.Ca_p);
Na_k = u(:, idx.Na_k);
K_k = u(:, idx.K_k);
Cl_k = u(:, idx.Cl_k);
HCO3_k = u(:, idx.HCO3_k);
Na_s = u(:, idx.Na_s);
K_s = u(:, idx.K_s);
HCO3_s = u(:, idx.HCO3_s);
w_k = u(:, idx.w_k);
I_k = u(:, idx.I_k);
Ca_k = u(:, idx.Ca_k);
h_k = u(:, idx.h_k);
s_k = u(:, idx.s_k);
m_k = u(:, idx.m_k); 
eet_k = u(:, idx.eet_k);
NO_k = u(:, idx.NO_k);
AA_k = u(:, idx.AA_k);

% SMC
Ca_i = u(:, idx.Ca_i);
s_i = u(:, idx.s_i);
v_i = u(:, idx.v_i);
w_i = u(:, idx.w_i);
I_i = u(:, idx.I_i);
NO_i = u(:, idx.NO_i);
E_b = u(:, idx.E_b);
E_6c = u(:, idx.E_6c);
cGMP_i = u(:, idx.cGMP_i);
H_i = u(:, idx.H_i);
AA_i = u(:, idx.AA_i);

% EC
Ca_j = u(:, idx.Ca_j);
s_j = u(:, idx.s_j);
v_j = u(:, idx.v_j);
I_j = u(:, idx.I_j);
eNOS_act_j = u(:, idx.eNOS_act_j);
NO_j = u(:, idx.NO_j);

% Wall Mechanics
Mp = u(:, idx.Mp);
AMp = u(:, idx.AMp);
AM = u(:, idx.AM);
R = u(:, idx.R);