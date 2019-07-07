%% Initial conditions for the NVU model

function u0 = initial_conditions(idx)

% Neuron
u0(idx.E_t) = 0;            % [-]
u0(idx.I_t) = 0;            % [-]
u0(idx.K_e) = 3.5;          % [mM]
u0(idx.Na_sa) = 9.37;       % [mM]
u0(idx.Na_d) = 9.42;        % [mM]
u0(idx.O2) = 0.02566;        % [uM]old value 0.0266
u0(idx.CBV) = 1.204;       % [-] old value 1.3167
u0(idx.HbR) = 0.7641;       % [-] old value 0.6665
u0(idx.Ca_n) = 0.1;         % [uM]
u0(idx.nNOS_act_n) = 0;%0.01056; % [uM]old value 0.318
u0(idx.NO_n) = 0.02425;      % [uM] old value 0.1671   

% Astrocyte
u0(idx.Na_k) = 18740;   % [uM] old value 18268 
u0(idx.K_k) = 92660;    % [uM]  old value 92708
u0(idx.HCO3_k) = 9085;  % [uM] old value 9131
u0(idx.Cl_k) = 8212;    % [uM] 0ld value 7733
u0(idx.Na_s) = 149200;  % [uM] old value 150255
u0(idx.K_s) = 2932;     % [uM] old vlaue 2837
u0(idx.HCO3_s) = 16980; % [uM] old value 16881
u0(idx.K_p) = 3039;   % [uM] old value 3045.1
u0(idx.w_k) = 8.26e-5; % [-] old value 1.703e-4
u0(idx.Ca_k) = 0.1435;  % [uM] old value 0.1612
u0(idx.s_k) = 480.8;    % [uM] 
u0(idx.h_k) = 0.4107;   % [-] old value 0.3828
u0(idx.I_k) = 0.048299; % [uM]
u0(idx.eet_k) = 0.4350; % [uM] old value 0.6123
u0(idx.m_k) = 0.513;   % [-] old value 0.5710
u0(idx.Ca_p) = 1853;  % [uM] old value 1746.4
u0(idx.NO_k) = 0.02234;  % [uM] old value 0.1106
u0(idx.v_k) = -88.79;      % [mV] 
u0(idx.AA_k) = 9.3; %[uM]

% SMC
u0(idx.Ca_i) = 0.2641; % 0.1;           [uM]
u0(idx.s_i) = 1.1686; % 0.1;            [uM]
u0(idx.v_i) = -34.7; % -60;             [mV]
u0(idx.w_i) = 0.2206; % 0.1;            [-]
u0(idx.I_i) = 0.275; % 0.1;             [uM]
u0(idx.NO_i) = 0.02047; % 0.05;          [uM] old value 0.0541
u0(idx.E_b) = 0.6372; % 1/3;            [-] old value 0.4077
u0(idx.E_6c) = 0.2606; % 1/3;           [-] old value 0.4396
u0(idx.cGMP_i) = 6.1; % 8;           [uM] old value 8.2826
u0(idx.H_i) = 0.096; % [uM] old value ).1271
u0(idx.AA_i) = 9.3; % [uM]

% EC
u0(idx.Ca_j) = 0.8339; % 0.1;           [uM]
u0(idx.s_j) = 0.6262; % 0.1;            [uM]
u0(idx.v_j) = -68.39; % -75;            [mV]
u0(idx.I_j) = 0.825; % 0.1;             [uM]
u0(idx.eNOS_act_j) = 0;%0.4451; % 0.7;     [uM]
u0(idx.NO_j) = 0.02051; % 0.05;          [uM] old value 0.0528

% Wall Mechanics
u0(idx.Mp) = 0.0842;                                                    % [-]
u0(idx.AMp) = 0.0622;                                                   % [-]
u0(idx.AM) = 0.2746; 
u0(idx.R) = 20;%22.44;

end