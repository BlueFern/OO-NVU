%% Indices for the variables in NVU model

function idx = all_indices()
% Neuron
idx.E_t = 1;             
idx.I_t = 2;           
idx.K_e = 3;           
idx.Na_sa = 4;            
idx.Na_d = 5;          
idx.O2 = 6;               
idx.CBV = 7;          
idx.HbR = 8;              
idx.Ca_n = 9;             
idx.nNOS_act_n = 10;  
idx.NO_n = 11;            

% Astrocyte
idx.v_k = 12;
idx.K_p = 13;
idx.Ca_p = 14;
idx.Na_k = 15;
idx.K_k = 16;
idx.Cl_k = 17;
idx.HCO3_k = 18;
idx.Na_s = 19;
idx.K_s = 20;
idx.HCO3_s = 21;
idx.w_k = 22;
idx.I_k = 23;
idx.Ca_k = 24;
idx.h_k = 25;
idx.s_k = 26;
idx.m_k = 27;
idx.eet_k = 28;
idx.NO_k = 29;
idx.AA_k = 30;

% SMC
idx.Ca_i = 31;
idx.s_i = 32;
idx.v_i = 33;
idx.w_i = 34;
idx.I_i = 35;
idx.NO_i = 36;
idx.E_b = 37;
idx.E_6c = 38;
idx.cGMP_i = 39;
idx.H_i = 40;
idx.AA_i = 41;

% EC
idx.Ca_j = 42;
idx.s_j = 43;
idx.v_j = 44;
idx.I_j = 45;
idx.eNOS_act_j = 46; 
idx.NO_j = 47;

% Wall Mechanics
idx.Mp = 48;
idx.AMp = 49;
idx.AM = 50;
idx.R = 51;

end