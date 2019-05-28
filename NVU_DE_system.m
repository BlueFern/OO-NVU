%% ODE system for NVU Version 3

function du = NVU_DE_system(t, u, idx, p, input_data) % input_data is used in all_algebraic_variables, don't remove

% Initialise size of du
du = zeros(p.num_variables,1);

% Set variable names 
E_t =   u(idx.E_t);        
I_t =   u(idx.I_t);            
K_e =   u(idx.K_e);           
Na_sa = u(idx.Na_sa);         
Na_d =  u(idx.Na_d);          
O2 =    u(idx.O2);             
CBV =   u(idx.CBV);           
HbR =   u(idx.HbR);           
Ca_n =  u(idx.Ca_n);          
nNOS_act_n = u(idx.nNOS_act_n); 
NO_n =  u(idx.NO_n);  
v_k = u(idx.v_k);
K_p = u(idx.K_p);
Ca_p = u(idx.Ca_p);
Na_k = u(idx.Na_k);
K_k = u(idx.K_k);
Cl_k = u(idx.Cl_k);
HCO3_k = u(idx.HCO3_k);
Na_s = u(idx.Na_s);
K_s = u(idx.K_s);
HCO3_s = u(idx.HCO3_s);
w_k = u(idx.w_k);
I_k = u(idx.I_k);
Ca_k = u(idx.Ca_k);
h_k = u(idx.h_k);
s_k = u(idx.s_k);
m_k = u(idx.m_k); 
eet_k = u(idx.eet_k);
NO_k = u(idx.NO_k);
AA_k = u(idx.AA_k);
Ca_i = u(idx.Ca_i);
s_i = u(idx.s_i);
v_i = u(idx.v_i);
w_i = u(idx.w_i);
I_i = u(idx.I_i);
NO_i = u(idx.NO_i);
E_b = u(idx.E_b);
E_6c = u(idx.E_6c);
cGMP_i = u(idx.cGMP_i);
H_i = u(idx.H_i);
AA_i = u(idx.AA_i);
Ca_j = u(idx.Ca_j);
s_j = u(idx.s_j);
v_j = u(idx.v_j);
I_j = u(idx.I_j);
eNOS_act_j = u(idx.eNOS_act_j);
NO_j = u(idx.NO_j);
Mp = u(idx.Mp);
AMp = u(idx.AMp);
AM = u(idx.AM);
R = u(idx.R);

% Load all fluxes and algebraic variables
all_algebraic_variables()

%% Set the RHS of the ODEs

du(idx.E_t) =  1 ./ p.tau_e .* (-E_t + (k_e - p.r_e .* E_t) .* s_arg_e);            % [ms^-1] Change in excitatory cells firing      
du(idx.I_t) =  1 ./ p.tau_i .* (-I_t + (k_i - p.r_i .* I_t) .* s_arg_i);            % [ms^-1] Change in inhibitory cells firing     

if p.gamma_switch == 0
    du(idx.K_e, :) = - p.beta_K_e * (K_e - p.K_eBase) + p.alpha_K_e * p.beta_K_e .* (abs(E_t - I_t) ./ p.EI_relative);  % [mM ms^-1]  Change in extracellular potassium concentration
else
    du(idx.K_e, :) = - p.beta_K_e * (K_e - p.K_eBase) + Gamma(t,p,E_t,I_t);
end

du(idx.Na_sa) = - p.beta_Na_sa * (Na_sa - p.Na_saBase) + p.alpha_Na_sa * p.beta_Na_sa .* (abs(E_t - I_t) ./ p.EI_relative);  % [mM ms^-1]   Change in soma/axon sodium concentration                  
du(idx.Na_d) =  - p.beta_Na_d * (Na_d - p.Na_dBase) +  p.alpha_Na_d * p.beta_Na_d .* (abs(E_t - I_t) ./ p.EI_relative); % [mM ms^-1]   Change in dendrite sodium concentration
du(idx.O2) = J_O2_vascular - J_O2_background - J_O2_pump;  % [uM ms^-1]  Change in oxygen concentration             
du(idx.CBV) = 1/(p.tau_MTT + p.tau_TAT) .* ( CBF/p.CBF_init  - CBV.^(1/p.d) ); % [ms^-1] Change in normalized cerebral blood volume         
du(idx.HbR) = 1/p.tau_MTT * ( CMRO2./CMRO2_init - HbR./CBV .* f_out );  % [ms^-1]  Change in normalized deoxyhemoglobin           
du(idx.Ca_n) = (I_Ca_tot / (2 * p.Farad * p.V_spine) - (p.k_ex * (Ca_n - p.Ca_rest))) / (1 + p.lambda_buf); % [uM ms^-1] Change in calcium concentration          
du(idx.nNOS_act_n) = p.V_maxNOS * CaM ./ (p.K_actNOS + CaM) - p.mu2_n * nNOS_act_n;         % [uM ms^-1]  Change in active NO concentration 
du(idx.NO_n) = p_NO_n - c_NO_n + d_NO_n;        % [uM ms^-1] Change in NO concentration         

du(idx.K_p) = J_BK_k ./ (p.VR_pa) + J_KIR_i ./ p.VR_ps - p.R_decay * (K_p - p.K_p_min) ;
du(idx.Ca_p) = J_TRPV_k./ p.VR_pa + J_VOCC_i ./ p.VR_ps - p.Ca_decay_k .* (Ca_p - p.Capmin_k); % calcium concentration in PVS
du(idx.v_k) = p.gamma_i .* ( -J_BK_k - J_K_k - J_Cl_k - J_NBC_k - J_Na_k - J_NaK_k - 2*J_TRPV_k);
du(idx.Na_k) = -J_Na_k - 3*J_NaK_k + J_NKCC1_k + J_NBC_k;
du(idx.K_k) = -J_K_k + 2*J_NaK_k + J_NKCC1_k + J_KCC1_k - J_BK_k;
du(idx.HCO3_k) = 2*J_NBC_k;
du(idx.Cl_k) = du(idx.Na_k) + du(idx.K_k) - du(idx.HCO3_k) + p.z_Ca * du(idx.Ca_k);

du(idx.K_s) = 1./p.VR_sa * (J_K_k - 2 * J_NaK_k - J_NKCC1_k - J_KCC1_k) + J_K_NEtoSC;
du(idx.Na_s) = 1./p.VR_sa * (J_Na_k + 3 * J_NaK_k - J_NKCC1_k - J_NBC_k) + J_Na_NEtoSC;
du(idx.HCO3_s) = 1./p.VR_sa * (-2*J_NBC_k);
du(idx.w_k) = phi_w .* (w_inf - w_k);
du(idx.I_k) =  p.r_h * G - p.k_deg * I_k;
du(idx.Ca_k) = B_cyt .* (J_IP3 - J_pump + J_ER_leak - J_TRPV_k/p.r_buff);    
du(idx.h_k) = p.k_on * (p.K_inh - (Ca_k + p.K_inh) .* h_k);
du(idx.s_k) = -(B_cyt .* (J_IP3 - J_pump + J_ER_leak)) ./ (p.VR_ER_cyt);
du(idx.m_k) =  p.trpv_switch .* ((minf_k - m_k) ./ p.t_TRPV_k) ; % TRPV open probability, the trpv_switch can switch the turn the channel on or off.
du(idx.eet_k) = p.V_eet * max(Ca_k - p.Ca_k_min, 0) - p.k_eet * eet_k;
du(idx.NO_k) = p_NO_k - c_NO_k + d_NO_k;
du(idx.AA_k) = (p.AA_m * p.AA_max)./(p.AA_m + (Ca_k - p.Ca0)).^2 .* du(idx.Ca_k) + (AA_i - AA_k)/p.tau_AA; 

du(idx.Ca_i) = J_IP3_i - J_SR_uptake_i - J_extrusion_i + J_SR_leak_i - J_VOCC_i + J_CICR_i + J_NaCa_i - 0.1*J_stretch_i + J_Ca_coup_i;
du(idx.s_i) = J_SR_uptake_i - J_CICR_i - J_SR_leak_i;
du(idx.v_i) = p.gamma_i * ( -J_NaK_i - J_Cl_i - 2*J_VOCC_i - J_NaCa_i - J_K_i -J_stretch_i - J_KIR_i) + V_coup_i;
du(idx.w_i) = p.lambda_i * (K_act_i - w_i); 
du(idx.I_i) = J_IP3_coup_i - J_degrad_i;
du(idx.NO_i) = p_NO_i - c_NO_i + d_NO_i;
du(idx.E_b) = -p.k1 * E_b .* NO_i + p.k_1 * E_6c + k4 .* E_5c;     
du(idx.E_6c) = p.k1 * E_b .* NO_i - (p.k_1 + p.k2) * E_6c - p.k3 * E_6c .* NO_i;
du(idx.cGMP_i) = p.V_max_sGC * E_5c - V_max_pde .* cGMP_i ./ (p.K_m_pde + cGMP_i); 
du(idx.H_i) = f_NO .* p.V_a .* AA_i ./(p.K_a + AA_i) +  p.V_f .* AA_i ./ (p.K_f + AA_i) - p.lambda_h * H_i; % 20-HETE increases with AA (but slower rate)
du(idx.AA_i) = (AA_k - AA_i)/p.tau_AA; % AA diffuses to SMC from the astrocyte

du(idx.Ca_j) = J_IP3_j - J_ER_uptake_j + J_CICR_j - J_extrusion_j + J_ER_leak_j + J_cation_j + p.J_0_j - J_stretch_j - J_Ca_coup_i;
du(idx.s_j) = J_ER_uptake_j - J_CICR_j - J_ER_leak_j; 
du(idx.v_j) = -1/p.C_m_j * (J_K_j + J_R_j) - V_coup_i;
du(idx.I_j) = p.J_PLC - J_degrad_j - J_IP3_coup_i;
du(idx.eNOS_act_j) = p.gam_eNOS .* Act_eNOS_Ca .* p.NOswitch_EC_CA  + (1 - p.gam_eNOS) .* Act_eNOS_wss .* p.NOswitch_EC_WSS - p.mu2_j * eNOS_act_j; 
du(idx.NO_j) = p_NO_j - c_NO_j + d_NO_j; 

du(idx.Mp) = p.wallMech * ( p.K_4 * AMp + K_1 .* M - (K_2 + p.K_3) .* Mp );
du(idx.AMp) = p.wallMech  * ( p.K_3 * Mp + K_6 .* AM - (p.K_4 + K_5) .* AMp );
du(idx.AM) = p.wallMech  * ( K_5 .* AMp - (p.K_7 + K_6) .* AM );
du(idx.R) =  p.R_init / eta * ( R * p.trans_p ./ h - E .* (R - R_0) ./ R_0);

end