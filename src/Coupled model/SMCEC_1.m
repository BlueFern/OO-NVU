classdef SMCEC_1 < handle
    % The 'SMCEC' code contains the following sections of the model:
    %   The Smooth Muscle cell and the Endothelial Cell
    %   Please refer to the relevant sections in the documentation for 
    %   full information on the equations and variable names.
    
    % SMC and EC of cell 1 of coupled model!
    
    %%%%% Both SMC and EC have coupling terms!!!
    %%%%% Coupling is set to 0 by default
    
    properties
        params
        u0
        index
        n_out
        idx_out
        enabled
    end
    methods
        function self = SMCEC_1(varargin)
            self.params = parse_inputs(varargin{:});
            self.index = indices();
            self.u0 = initial_conditions(self.index);
            self.enabled = true(size(self.u0));
            [self.idx_out, self.n_out] = output_indices();
        end
        function [du, varargout] = rhs(self, t, u, R_1, h_1, K_p_1, Ca_i_2, I_i_2, v_i_2, Ca_j_2, I_j_2, v_j_2)
            % Initalise inputs and parameters
            p = self.params;
            idx = self.index;
            
            Ca_i_1 = u(idx.Ca_i_1, :);
            s_i_1 = u(idx.s_i_1, :);
            v_i_1 = u(idx.v_i_1, :);
            w_i_1 = u(idx.w_i_1, :);
            I_i_1 = u(idx.I_i_1, :);
            Ca_j_1 = u(idx.Ca_j_1, :);
            s_j_1 = u(idx.s_j_1, :);
            v_j_1 = u(idx.v_j_1, :);
            I_j_1 = u(idx.I_j_1, :);
            
            %% SMC fluxes
            J_IP3_i_1 = p.F_i * I_i_1.^2 ./ (p.K_r_i^2 + I_i_1.^2);
            J_SR_uptake_i_1 = p.B_i * Ca_i_1.^2 ./ (p.c_b_i^2 + Ca_i_1.^2);
            J_CICR_i_1 = p.C_i * s_i_1.^2 ./ (p.s_c_i^2 + s_i_1.^2) .* ...
                Ca_i_1.^4 ./ (p.c_c_i^4 + Ca_i_1.^4);
            J_extrusion_i_1 = p.D_i * Ca_i_1 .* (1 + (v_i_1 - p.v_d) / p.R_d_i);
            J_SR_leak_i_1 = p.L_i * s_i_1;
            J_VOCC_i_1 = p.G_Ca_i * (v_i_1 - p.v_Ca1_i) ./ ...
                (1 + exp(-(v_i_1 - p.v_Ca2_i) / p.R_Ca_i));
            J_NaCa_i_1 = p.G_NaCa_i * Ca_i_1 ./ (Ca_i_1 + p.c_NaCa_i) .* ...
                (v_i_1 - p.v_NaCa_i);
            J_stretch_i_1 = p.G_stretch ./ ...
                (1 + exp(-p.alpha_stretch*(p.delta_p*R_1./h_1 - p.sigma_0))) .* ...
                (v_i_1 - p.E_SAC);
            J_NaK_i_1 = p.F_NaK_i;
            J_Cl_i_1 = p.G_Cl_i * (v_i_1 - p.v_Cl_i);
            J_K_i_1 = p.G_K_i * w_i_1 .* (v_i_1 - p.v_K_i);
            
            J_KIR_i_1 = self.shared(t, u, K_p_1);
            
            J_degrad_i_1 = p.k_d_i * I_i_1;
            
            %% SMC coupling
            J_Ca_SMCcoup_i_1 = p.D_Ca_i * (Ca_i_2 - Ca_i_1);
            J_IP3_SMCcoup_i_1 = p.D_IP3_i * (I_i_2 - I_i_1);
            J_V_SMCcoup_i_1 = p.D_v_i * (v_i_2 - v_i_1);
            
            %% EC coupling
            J_Ca_ECcoup_j_1 = p.D_Ca_j * (Ca_j_2 - Ca_j_1);
            J_IP3_ECcoup_j_1 = p.D_IP3_j * (I_j_2 - I_j_1);
            J_V_ECcoup_j_1 = p.D_v_j * (v_j_2 - v_j_1);
            
            %% EC fluxes
            J_IP3_j_1 = p.F_j * I_j_1.^2 ./ (p.K_r_j^2 + I_j_1.^2);
            J_ER_uptake_j_1 = p.B_j * Ca_j_1.^2 ./ (p.c_b_j^2 + Ca_j_1.^2);
            J_CICR_j_1 = p.C_j * s_j_1.^2 ./ (p.s_c_j^2 + s_j_1.^2) .* ...
                Ca_j_1.^4 ./ (p.c_c_j^4 + Ca_j_1.^4);
            J_extrusion_j_1 = p.D_j * Ca_j_1;
            J_stretch_j_1 = p.G_stretch ./ ...
                (1 + exp(-p.alpha_stretch*(p.delta_p*R_1./h_1 - p.sigma_0))) .* ...
                (v_j_1 - p.E_SAC);
            
            J_ER_leak_j_1 = p.L_j * s_j_1;
            
            J_cation_j_1 = p.G_cat_j * (p.E_Ca_j - v_j_1) * 0.5 .* ...
                (1 + tanh((log10(Ca_j_1) - p.m_3_cat_j) / p.m_4_cat_j));
            
            J_BK_Ca_j_1 = 0.2 * (1 + tanh( ...
                ((log10(Ca_j_1) - p.c) .* (v_j_1 - p.bb_j) - p.a_1_j) ./ ...
                (p.m_3b_j * (v_j_1 + p.a_2_j*(log10(Ca_j_1) - p.c) - p.bb_j).^2 + ...
                p.m_4b_j)));
            J_SK_Ca_j_1 = 0.3 * (1 + tanh((log10(Ca_j_1) - p.m_3s_j) / p.m_4s_j));
            J_K_j_1 = p.G_tot_j * (v_j_1 - p.v_K_j) .* (J_BK_Ca_j_1 + J_SK_Ca_j_1);
            J_R_j_1 = p.G_R_j * (v_j_1 - p.v_rest_j);
            J_degrad_j_1 = p.k_d_j * I_j_1;
            
            %% Coupling between SMC and EC
            V_coup_i_1 = -p.G_coup * (v_i_1 - v_j_1);
            J_IP3_coup_i_1 = -p.P_IP3 * (I_i_1 - I_j_1);
            J_Ca_coup_i_1 = -p.P_Ca * (Ca_i_1 - Ca_j_1);
            
            K_act_i_1 = (Ca_i_1 + p.c_w_i).^2 ./ ...
                ((Ca_i_1 + p.c_w_i).^2 + ...
                p.beta_i * exp(-(v_i_1 - p.v_Ca3_i) / p.R_K_i));
            
            
            %% Differential Equations
            % Smooth muscle cell
            du(idx.Ca_i_1, :) = J_IP3_i_1 - J_SR_uptake_i_1 - J_extrusion_i_1 + ...
                J_SR_leak_i_1 - J_VOCC_i_1 + J_CICR_i_1 + J_NaCa_i_1 + ...
                0.1*J_stretch_i_1 + J_Ca_coup_i_1 + J_Ca_SMCcoup_i_1;
            
            du(idx.s_i_1, :) = J_SR_uptake_i_1 - J_CICR_i_1 - J_SR_leak_i_1;
            du(idx.v_i_1, :) = p.gamma_i * (...
                -J_NaK_i_1 - J_Cl_i_1 - 2*J_VOCC_i_1 - J_NaCa_i_1 - J_K_i_1 ...
                -J_stretch_i_1 - J_KIR_i_1) + V_coup_i_1 + J_V_SMCcoup_i_1;
            du(idx.w_i_1, :) = p.lambda_i * (K_act_i_1 - w_i_1);
            du(idx.I_i_1, :) = J_IP3_coup_i_1 - J_degrad_i_1 + J_IP3_SMCcoup_i_1;
            du(idx.K_i_1, :) = J_NaK_i_1 - J_KIR_i_1 - J_K_i_1;
            % Endothelial Cell
            du(idx.Ca_j_1, :) = J_IP3_j_1 - J_ER_uptake_j_1 + J_CICR_j_1 - ...
                J_extrusion_j_1 + J_ER_leak_j_1 + J_cation_j_1 + p.J_0_j_1 + ...
                J_stretch_j_1 - J_Ca_coup_i_1 + J_Ca_ECcoup_j_1;
            du(idx.s_j_1, :) = J_ER_uptake_j_1 - J_CICR_j_1 - J_ER_leak_j_1;
            du(idx.v_j_1, :) = -1/p.C_m_j * (J_K_j_1 + J_R_j_1) - V_coup_i_1 + J_V_ECcoup_j_1;
            du(idx.I_j_1, :) = p.J_PLC_1 - J_degrad_j_1 - J_IP3_coup_i_1 + J_IP3_ECcoup_j_1;
            
            du = bsxfun(@times, self.enabled, du);
            
            if nargout == 2
                Uout = zeros(self.n_out, size(u, 2));
                Uout(self.idx_out.V_coup_i_1, :) = V_coup_i_1;
                Uout(self.idx_out.J_Ca_coup_i_1, :) = J_Ca_coup_i_1;
                Uout(self.idx_out.J_IP3_coup_i_1, :) = J_IP3_coup_i_1;
                Uout(self.idx_out.J_stretch_i_1, :) = J_stretch_i_1;
                Uout(self.idx_out.J_IP3_i_1, :) = J_IP3_i_1;
                Uout(self.idx_out.J_SR_uptake_i_1, :) = J_SR_uptake_i_1;
                Uout(self.idx_out.J_CICR_i_1, :) = J_CICR_i_1;
                Uout(self.idx_out.J_extrusion_i_1, :) = J_extrusion_i_1;
                Uout(self.idx_out.J_SR_leak_i_1, :) = J_SR_leak_i_1;
                Uout(self.idx_out.J_VOCC_i_1, :) = J_VOCC_i_1;
                Uout(self.idx_out.J_NaCa_i_1, :) = J_NaCa_i_1;
                Uout(self.idx_out.J_NaK_i_1, :) = J_NaK_i_1;
                Uout(self.idx_out.J_Cl_i_1, :) = J_Cl_i_1;
                Uout(self.idx_out.J_K_i_1, :) = J_K_i_1;
                Uout(self.idx_out.J_KIR_i_1, :) = J_KIR_i_1;
                Uout(self.idx_out.K_act_i_1, :) = K_act_i_1;
                Uout(self.idx_out.J_degrad_i_1, :) = J_degrad_i_1;
                
                Uout(self.idx_out.V_coup_j_1, :) = -V_coup_i_1;
                Uout(self.idx_out.J_Ca_coup_j_1, :) = -J_Ca_coup_i_1;
                Uout(self.idx_out.J_IP3_coup_j_1, :) = -J_IP3_coup_i_1;
                Uout(self.idx_out.J_stretch_j_1, :) = J_stretch_j_1;
                Uout(self.idx_out.J_0_j_1, :) = p.J_0_j_1;
                Uout(self.idx_out.J_IP3_j_1, :) = J_IP3_j_1;
                Uout(self.idx_out.J_ER_uptake_j_1, :) = J_ER_uptake_j_1;
                Uout(self.idx_out.J_CICR_j_1, :) = J_CICR_j_1;
                Uout(self.idx_out.J_extrusion_j_1, :) = J_extrusion_j_1;
                Uout(self.idx_out.J_ER_leak_j_1, :) = J_ER_leak_j_1;
                Uout(self.idx_out.J_cation_j_1, :) = J_cation_j_1;
                Uout(self.idx_out.J_BK_Ca_j_1, :) = J_BK_Ca_j_1;
                Uout(self.idx_out.J_SK_Ca_j_1, :) = J_SK_Ca_j_1;
                Uout(self.idx_out.J_K_j_1, :) = J_K_j_1;
                Uout(self.idx_out.J_R_j_1, :) = J_R_j_1;
                Uout(self.idx_out.J_degrad_j_1, :) = J_degrad_j_1;
                
                Uout(self.idx_out.J_Ca_SMCcoup_i_1, :) = J_Ca_SMCcoup_i_1;
                Uout(self.idx_out.J_IP3_SMCcoup_i_1, :) = J_IP3_SMCcoup_i_1;
                Uout(self.idx_out.J_V_SMCcoup_i_1, :) = J_V_SMCcoup_i_1;
                
                Uout(self.idx_out.J_Ca_ECcoup_j_1, :) = J_Ca_ECcoup_j_1;
                Uout(self.idx_out.J_IP3_ECcoup_j_1, :) = J_IP3_ECcoup_j_1;
                Uout(self.idx_out.J_V_ECcoup_j_1, :) = J_V_ECcoup_j_1;
                varargout{1} = Uout; 
            end
        end
        function [J_KIR_i_1, Ca_i_1, I_i_1, v_i_1, Ca_j_1, I_j_1, v_j_1] = shared(self, ~, u, K_p_1)
            p = self.params;
            idx = self.index;
            v_i_1 = u(idx.v_i_1, :);
            Ca_i_1 = u(idx.Ca_i_1, :);
            I_i_1 = u(idx.I_i_1, :);
            v_j_1 = u(idx.v_j_1, :);
            Ca_j_1 = u(idx.Ca_j_1, :);
            I_j_1 = u(idx.I_j_1, :);
            v_KIR_i = p.z_1 * K_p_1 - p.z_2;
            g_KIR_i = exp(p.z_5 * v_i_1 + p.z_3 * K_p_1 - p.z_4);
            J_KIR_i_1 = p.F_KIR_i * g_KIR_i / p.gamma_i .* (v_i_1 - v_KIR_i);
        end
        function names = varnames(self)
            names = [fieldnames(self.index); fieldnames(self.idx_out)];
        end
        
        
    end
end

function idx = indices()
% Index of parameters needing inital conditions 
idx.Ca_i_1 = 1;
idx.s_i_1 = 2;
idx.v_i_1 = 3;
idx.w_i_1 = 4;
idx.I_i_1 = 5;
idx.K_i_1 = 6;
idx.Ca_j_1 = 7;
idx.s_j_1 = 8;
idx.v_j_1 = 9;
idx.I_j_1 = 10;
end

function [idx, n] = output_indices()
% Index of all other output parameters
idx.V_coup_i_1 = 1;
idx.J_Ca_coup_i_1 = 2;
idx.J_IP3_coup_i_1 = 3;
idx.J_stretch_i_1 = 4;
idx.J_IP3_i_1 = 5;
idx.J_SR_uptake_i_1 = 6;
idx.J_CICR_i_1 = 7;
idx.J_extrusion_i_1 = 8;
idx.J_SR_leak_i_1 = 9;
idx.J_VOCC_i_1 = 10;
idx.J_NaCa_i_1 = 11;
idx.J_NaK_i_1 = 12;
idx.J_Cl_i_1 = 13;
idx.J_K_i_1 = 14;
idx.J_KIR_i_1 = 15;
idx.K_act_i_1 = 16;
idx.J_degrad_i_1 = 17;

idx.V_coup_j_1 = 18;
idx.J_Ca_coup_j_1 = 19;
idx.J_IP3_coup_j_1 = 20;
idx.J_stretch_j_1 = 21;
idx.J_0_j_1 = 22;
idx.J_IP3_j_1 = 23;
idx.J_ER_uptake_j_1 = 24;
idx.J_CICR_j_1 = 25;
idx.J_extrusion_j_1 = 26;
idx.J_ER_leak_j_1 = 27;
idx.J_cation_j_1 = 28;
idx.J_BK_Ca_j_1 = 29;
idx.J_SK_Ca_j_1 = 30;
idx.J_K_j_1 = 31;
idx.J_R_j_1 = 32;
idx.J_degrad_j_1 = 33;

idx.J_Ca_SMCcoup_i_1 = 34;
idx.J_IP3_SMCcoup_i_1 = 35;
idx.J_V_SMCcoup_i_1 = 36;

idx.J_Ca_ECcoup_j_1 = 37;
idx.J_IP3_ECcoup_j_1 = 38;
idx.J_V_ECcoup_j_1 = 39;

n = numel(fieldnames(idx));
end

function params = parse_inputs(varargin)
parser = inputParser();
% Coupling constants for SMC
parser.addParameter('D_Ca_i', 0); % s^-1                          %%%%%%%%%
parser.addParameter('D_IP3_i', 0); % s^-1                         %%%%%%%%%
parser.addParameter('D_v_i', 0); % s^-1                           %%%%%%%%%

% Coupling constants for EC
parser.addParameter('D_Ca_j', 0); % s^-1                          %%%%%%%%%
parser.addParameter('D_IP3_j', 0); % s^-1                         %%%%%%%%%
parser.addParameter('D_v_j', 0); % s^-1                           %%%%%%%%%

% Smooth Muscle Cell ODE Constants
parser.addParameter('gamma_i', 1970); %mV uM^-1
parser.addParameter('lambda_i', 45); % s^-1

% Endothelial Cell ODE Constants
parser.addParameter('C_m_j', 25.8); %pF
parser.addParameter('J_PLC_1', 0.18); % uMs^-1                  %%%%%%%%%
parser.addParameter('J_0_j_1', 0.029); %constant Ca influx (EC)

% Smooth Muscle Cell Flux Constants
parser.addParameter('F_i', 0.23); %uM s^-1
parser.addParameter('K_r_i', 1); %uM

parser.addParameter('B_i', 2.025); %uM s^-1
parser.addParameter('c_b_i', 1.0); %uM

parser.addParameter('C_i', 55); % uM s^-1
parser.addParameter('s_c_i', 2.0); %uM
parser.addParameter('c_c_i', 0.9); %uM

parser.addParameter('D_i', 0.24); %s^-1
parser.addParameter('v_d', -100); %mV
parser.addParameter('R_d_i', 250); %mV

parser.addParameter('L_i', 0.025); %s^-1

parser.addParameter('G_Ca_i', 1.29e-3); %uM mV^-1 s^-1
parser.addParameter('v_Ca1_i', 100); %mV
parser.addParameter('v_Ca2_i', -24); %mV
parser.addParameter('R_Ca_i', 8.5); %mV

parser.addParameter('G_NaCa_i', 3.16e-3); %uM mV^-1 s^-1
parser.addParameter('c_NaCa_i', 0.5); %uM
parser.addParameter('v_NaCa_i', -30); %mV

parser.addParameter('G_stretch', 6.1e-3); % uM mV^-1 s^-1   (Also EC parameter)
parser.addParameter('alpha_stretch', 7.4e-3); % mmHg^-1     (Also EC parameter)
parser.addParameter('delta_p', 30); % mmHg                  (Also EC parameter)
parser.addParameter('sigma_0', 500); % mmHg                 (Also EC parameter)
parser.addParameter('E_SAC', -18); % mV                     (Also EC parameter)

parser.addParameter('F_NaK_i', 4.32e-2); %uM s^-1

parser.addParameter('G_Cl_i', 1.34e-3); %uM mV^-1 s^-1
parser.addParameter('v_Cl_i', -25); %mV

parser.addParameter('G_K_i', 4.46e-3); %uM mV^-1 s^-1
parser.addParameter('v_K_i', -94); %mV

parser.addParameter('F_KIR_i', 7.5e2); % [-]
parser.addParameter('k_d_i', 0.1); % s^-1

% Endothelial Cell Flux Constants
parser.addParameter('F_j', 0.23); %uM s^-1
parser.addParameter('K_r_j', 1); %uM
parser.addParameter('B_j', 0.5); %uM s^-1
parser.addParameter('c_b_j', 1); %uM
parser.addParameter('C_j', 5); %uM s^-1
parser.addParameter('s_c_j', 2); %uM
parser.addParameter('c_c_j', 0.9); %uM
parser.addParameter('D_j', 0.24);% s^-1

% (G_stretch, alpha_stretch, delta_p, sigma0, E_SAC are included above in 
%  SMC flux Constants)

parser.addParameter('L_j', 0.025); %s^-1 

parser.addParameter('G_cat_j', 6.6e-4); %uM mV^-1 s^-1
parser.addParameter('E_Ca_j', 50); %mV
parser.addParameter('m_3_cat_j', -0.18); %uM
parser.addParameter('m_4_cat_j', 0.37); %uM

parser.addParameter('G_tot_j', 6927); %p mho
parser.addParameter('v_K_j', -80); %m mho

parser.addParameter('c', -0.4); %uM
parser.addParameter('bb_j', -80.8); %mV
parser.addParameter('a_1_j', 53.3); %uM mV
parser.addParameter('a_2_j', 53.3); % mV uM^-1
parser.addParameter('m_3b_j', 1.32e-3); %uM mV^-1
parser.addParameter('m_4b_j', 0.3); %uM mV
parser.addParameter('m_3s_j', -0.28); %uM
parser.addParameter('m_4s_j', 0.389); %uM

parser.addParameter('G_R_j', 955); %p omh
parser.addParameter('v_rest_j', -31.1); %mV
parser.addParameter('k_d_j', 0.1); %s^-1

parser.addParameter('P_Ca', 0.05); %s^-1
parser.addParameter('P_IP3', 0.05); %s^-1
parser.addParameter('G_coup', 0.5); %s^-1

% Additional Equations Constants
parser.addParameter('c_w_i', 0); %uM
parser.addParameter('beta_i', 0.13); %uM^2
parser.addParameter('v_Ca3_i', -27); %mV
parser.addParameter('R_K_i', 12); %mV
parser.addParameter('z_1', 4.5e-3); %mV
parser.addParameter('z_2', 112); %mV
parser.addParameter('z_3', 4.2e-4); %uM mV^-1 s^-1
parser.addParameter('z_4', 12.6); %uM mV^-1 s^-1
parser.addParameter('z_5', -7.4e-2); %uM mV^-1 s^-1

parser.parse(varargin{:})
params = parser.Results;
end

function u0 = initial_conditions(idx)
u0 = zeros(length(fieldnames(idx)), 1);
% Inital estimations of parameters from experimental data
u0(idx.Ca_i_1) = 0.1;
u0(idx.s_i_1) = 0.1;
u0(idx.v_i_1) = -60;
u0(idx.w_i_1) = 0.1;
u0(idx.I_i_1) = 0.1;

u0(idx.K_i_1) = 100e3;

u0(idx.Ca_j_1) = 0.1;
u0(idx.s_j_1) = 0.1;
u0(idx.v_j_1) = -75;
u0(idx.I_j_1) = 0.1;
end


