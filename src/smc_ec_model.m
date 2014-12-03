function [f_rhs, idx, p] = smc_ec_model(varargin)
idx = indices();
p = parse_inputs(varargin{:});

f_rhs = @rhs;

    function du = rhs(~, u, R, h, K_p)
        Ca_i = u(idx.Ca_i, :);
        s_i = u(idx.s_i, :);
        v_i = u(idx.v_i, :);
        w_i = u(idx.w_i, :);
        I_i = u(idx.I_i, :);
        % K_i = u(idx.K_i, :); not currently used anywhere
        Ca_j = u(idx.Ca_j, :);
        s_j = u(idx.s_j, :);
        v_j = u(idx.v_j, :);
        I_j = u(idx.I_j, :);
        
        % SMC fluxes
        J_IP3_i = p.F_i * I_i.^2 ./ (p.K_r_i^2 + I_i.^2);
        J_SR_uptake_i = p.B_i * Ca_i.^2 ./ (p.c_b_i^2 + Ca_i.^2);
        J_CICR_i = p.C_i * s_i.^2 ./ (p.s_c_i^2 + s_i.^2) .* ...
            Ca_i.^4 ./ (p.c_c_i^4 + Ca_i.^4);
        J_extrusion_i = p.D_i * Ca_i .* (1 + (v_i - p.v_d) / p.R_d_i);
        J_SR_leak_i = p.L_i * s_i;
        J_VOCC_i = p.G_Ca_i * (v_i - p.v_Ca1_i) ./ ...
            (1 + exp(-(v_i - v_Ca2_i) / p.R_Ca_i));
        J_NaCa_i = p.G_NaCa_i * Ca_i ./ (Ca_i + p.c_NaCa_i) .* ...
            (v_i - p.v_NaCa_i);
        J_stretch_i = p.G_stretch ./ ...
            (1 + exp(-p.alpha_stretch*(p.delta_p*R./h - p.sigma_0))) .* ...
            (v_i - p.E_SAC);
        J_NaK_i = p.F_NaK_i;
        J_Cl_i = p.G_Cl_i * (v_i - p.v_Cl_i);
        J_K_i = p.G_K_u * w_i .* (v_i - p.v_K_i);
        
        
        v_KIR_i = p.z_1 * K_p - p.z_2;
        g_KIR_i = exp(p.z_5 * v_i + p.z_3 * K_p - p.z_4);
        
        J_KIR_i = p.F_KIR_i * g_KIR_i / p.gamma_i * (v_i - v_KIR_i);
        J_degrad_i = p.k_d_i * I_i;
        
        % EC fluxes
        J_IP3_j = p.F_j * I_j.^2 ./ (p.K_r_j^2 + I_j.^2);
        J_ER_uptake_j = p.B_j * Ca_j.^2 ./ (p.c_b_j^2 + Ca_j.^2);
        J_CICR_j = p.C_j * s_j.^2 ./ (p.s_c_j^2 + s_j.^2) .* ...
            Ca_j.^4 ./ (p.c_c_j^4 + Ca_j.^4);
        J_extrusion_j = p.D_j * Ca_j;
        J_stretch_j = p.G_stretch ./ ...
            (1 + exp(-p.alpha_stretch*(p.delta_p*R./h - p.sigma_0))) .* ...
            (v_j - v.E_SAC);
        
        J_ER_leak_j = p.L_j * s_j;
        
        J_cation_j = p.G_cat_j * (p.E_Ca_j - v_j) * 0.5 .* ...
            (1 + tanh((log10(Ca_j) - p.m_3_cat_j) / p.m_4_cat_j));
        
        J_BK_Ca_j = 0.2 * (1 + tanh( ...
            ((log10(Ca_j) - p.c) .* (v_j - p.b_j) - p.a_1_j) ./ ...
            (p.m_3b_j * (v_j + p.a_2_j*(log10(Ca_j) - p.c) - p.b_j)^2 + ...
            p.m_4b_j)));
        J_SK_Ca_j = 0.3 * (1 + tanh((log10(Ca_j) - p.m_3s_j) / p.m_4s_j));
        J_K_j = p.G_tot_j * (v_j - v_K_j) .* (J_BK_Ca_j + J_SK_Ca_j);
        J_R_j = p.G_R_j * (v_j - p.v_rest_j);
        J_degrad_j = p.k_d_j * I_j;
        
        % Coupling
        V_coup_i = -p.G_coup * (v_i - v_j);
        J_IP3_coup_i = -p.P_IP3 * (I_i - I_j);
        J_Ca_coup_i = -p.P_Ca * (Ca_i - Ca_j);
        
        K_act_i = (Ca_i + p.c_w_i).^2 ./ ...
            ((Ca_i + p.c_w_i).^2 + ...
            p.beta_i * exp(-(v_i - p.v_Ca3_i) / p.R_K_i));
        
        
        %ODE RHSes
        du(idx.Ca_i, :) = J_IP3_i - J_SR_uptake_i - J_extrusion_i + ...
            J_SR_leak_i - J_VOCC_i + J_NaCa_i + 0.1*J_stretch_i + ...
            J_Ca_coup_i;
        du(idx.s_i, :) = J_SR_uptake_i - J_CICR_i - J_RS_leak_i;
        du(idx.v_i, :) = p.gamma_i * (...
            -J_NaKa_i - J_Cl_i - 2*J_VOCC_i - J_NaCa_i - J_K_i ...
            -J_stretch_i - J_KIR_i) + V_coup_i;
        du(idx.w_i, :) = p.lambda_i * (K_act_i - w_i);
        du(idx.I_i, :) = J_IP3_coup_i - J_degrad_i;
        du(idx.K_i, :) = J_NaK_i - J_KIR_i - J_K_i;
        du(idx.Ca_j, :) = J_IP3_j - J_ER_uptake_j + J_CICR_j - ...
            J_extrusion_j + J_ER_leak_j + J_cation_j + p.J_0_j + ...
            J_stretch_j - J_Ca_coup_i;
        du(idx.s_j, :) = J_ER_uptake_j - J_CICR_j - J_ER_leak_j;
        du(idx.v_j, :) = -1/p.C_m_j * (J_K_j + J_R_j) - V_coup_i;
        du(idx.I_j, :) = J_PLC - J_degrad_j - J_IP3_coup_i;
    end
end
function idx = indices()
idx.Ca_i = 1;
idx.s_i = 2;
idx.v_i = 3;
idx.w_i = 4;
idx.I_i = 5;
idx.K_i = 6;
idx.Ca_j = 7;
idx.s_j = 8;
idx.v_j = 9;
idx.I_j = 10;
end

function params = parse_inputs(varargin)
parser = inputParser();
parser.addParameter('gamma_i', 1970); %mV uM^-1
parser.addParameter('lambda_i', 45); % s^-1
parser.addParameter('C_m_j', 25.8); %pF
parser.addParameter('F_i', 0.23); %uM s^-1
parser.addParameter('K_r_i', 1); %uM
parser.addParameter('B_i', 2.025); %uM s^-1
parser.addParameter('c_bi', 1.0); %uM
parser.addParameter('C_i', 55); % uM s^-1
parser.addParameter('s_ci', 2.0); %uM
parser.addParameter('c_ci', 0.9); %uM
parser.addParameter('D_i', 0.24); %s^-1
parser.addParameter('v_d', -100); %mV
parser.addParameter('R_di', 250); %mV
parser.addParameter('L_i', 0.025); %s^-1
parser.addParameter('G_Ca_i', 1.29e-3); %uM mV^-1 s^-1
parser.addParameter('v_Ca1_i', 100); %mV
parser.addParameter('v_Ca2_i', -24); %mV
parser.addParameter('R_Ca_i', 8.5); %mV
parser.addParameter('G_NaCa_i', 3.16e-3); %uM mV^-1 s^-1
parser.addParameter('c_NaCa_i', 0.5); %uM
parser.addParameter('v_NaCa_i', -30); %mV

parser.addParameter('G_stretch', 61e-3); % uM mV^-1 s^-1
parser.addParameter('alpha_stretch', 74e-3); % mmHg
parser.addParameter('delta_p', 30); % mmHg
parser.addParameter('sigma_0', 500); % mmHg
parser.addParameter('E_SAC', -18); % mV

parser.addParameter('F_NaK', 4.23e-3); %uM s^-1
parser.addParameter('G_Cl_i', 1.34e-3); %uM mV^-1 s^-1
parser.addParameter('v_Cl_i', -25); %mV

parser.addParameter('G_K_i', 4.46e-3); %uM mV^-1 s^-1
parser.addParameter('v_K_i', -94); %mV

parser.addParameter('F_KIR_i', 7.5e2);
parser.addParameter('k_d_i', 0.1); % s^-1
% Endothelial

parser.addParameter('F_j', 0.23); %uM s^-1
parser.addParameter('K_r_j', 1); %uM
parser.addParameter('B_j', 0.5); %uM s^-1
parser.addParameter('c_b_j', 1); %uM
parser.addParameter('C_j', 5); %uM s^-1
parser.addParameter('s_c_j', 2); %uM
parser.addParameter('c_c_j', 0.9); %uM
parser.addParameter('D_j', 0.24);% s^-1
parser.addParameter('L_j', 0.025); %s^-1

parser.addParameter('G_cat_j', 6.6e-4); %uM mV^-1 s^-1
parser.addParameter('E_Ca_j', 50); %mV
parser.addParameter('m_3_cat_j', -6.18); %uM
parser.addParameter('m_4_cat_j', 0.37); %uM

parser.addParameter('G_tot_j', 6927); %pS
parser.addParameter('v_K_j', -80); %mV

parser.addParameter('c', -6.4); %uM
parser.addParameter('b_j', -80.8); %mV
parser.addParameter('a_1_j', 53.3); %uM mV
parser.addParameter('a_2_j', 53.3); %..
parser.addParameter('m_3b_j', 1.32e-3); %uM mV^-1
parser.addParameter('m_4b_j', 0.3); %uM mV
parser.addParameter('m_3s_j', -0.28); %uM
parser.addParameter('m_4s_j', 0.389); %uM

parser.addParameter('G_R_j', 955); %pS
parser.addParameter('v_rest_j', -31.1); %mV
parser.addParameter('k_d_j', 0.1); %s^-1

parser.addParameter('P_Ca', 0.05); %s^-1
parser.addParameter('P_IP3', 0.05); %s^-1
parser.addParameter('G_coup', 5); %s^-1

parser.addParameter('c_w_i', 0); %uM
parser.addParameter('beta_i', 0.13); %uM^2
parser.addParameter('v_Ca3_i', -27); %mV
parser.addParameter('R_K_i', 12); %mV
parser.addParameter('z_1', 4.5e3); %mV
parser.addParameter('z_2', 112); %..
parser.addParameter('z_3', 4.2e2); %uM mV^-1 s^-1
parser.addParameter('z_4', 12.6); %uM mV^-1 s^-1
parser.addParameter('z_5', -7.4e-2); %uM mV^-1 s^-1

parser.addParameter('J_0_j', 0.029); %constant Ca influx (EC)

parser.addParameter('J_PLC', 0.4); 


parser.parse(varargin{:})
params = parser.Results;
end
