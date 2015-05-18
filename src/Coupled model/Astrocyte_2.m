classdef Astrocyte_2 < handle
    % The 'Astrocyte' code contains the following sections of the model:
    %   The Neuron, Synaptic cleft, the Astrocyte and the Perivascular Space
    %   Currently there is no content under the Neuron sub-section
    %   Please refer to the relevient sections in the documentation for 
    %   full information on the equations and variable names.
    %   Test commit
    properties
        params
        u0
        index
        n_out
        idx_out
        enabled
    end
    methods
        function self = Astrocyte_2(varargin)
            self.params = parse_inputs(varargin{:});
            self.index = indices();
            self.u0 = initial_conditions(self.index);
            self.enabled = true(size(self.u0));
            [self.idx_out, self.n_out] = output_indices();
        end
        function [du, varargout] = rhs(self, t, u, J_KIR_i_2)
            % Initalise inputs and parameters
            t = t(:).';
            p = self.params;
            idx = self.index;
            
            R_k_2 = u(idx.R_k_2, :);
            K_p_2 = u(idx.K_p_2, :);
            
            N_Na_k_2 = u(idx.N_Na_k_2, :);
            N_K_k_2 = u(idx.N_K_k_2, :);
            N_Cl_k_2 = u(idx.N_Cl_k_2, :);
            N_HCO3_k_2 = u(idx.N_HCO3_k_2, :);
            
            N_Na_s_2 = u(idx.N_Na_s_2, :);
            N_K_s_2 = u(idx.N_K_s_2, :);
            N_HCO3_s_2 = u(idx.N_HCO3_s_2, :);
                   
            w_k_2 = u(idx.w_k_2, :);

            du = zeros(size(u));
            
            %% Synaptic Cleft: 
            % Electroneutrality condition 
            N_Cl_s_2 = N_Na_s_2 + N_K_s_2 - N_HCO3_s_2;
            
            % Volume-surface ratio
            R_s_2 = p.R_tot - R_k_2; 
            
            % Scale concentrations to get actual concentrations,
            K_s_2 = N_K_s_2 ./ R_s_2;
            Na_s_2 = N_Na_s_2 ./ R_s_2;
            Cl_s_2 = N_Cl_s_2 ./ R_s_2;
            HCO3_s_2 = N_HCO3_s_2 ./ R_s_2;
            Na_k_2 = N_Na_k_2 ./ R_k_2;
            K_k_2 = N_K_k_2 ./ R_k_2;
            Cl_k_2 = N_Cl_k_2 ./ R_k_2;
            HCO3_k_2 = N_HCO3_k_2 ./ R_k_2;
            
            %% Astrocyte 
            % Volume-surface Ratio; Scaling ODE
            du(idx.R_k_2, :) = p.L_p * ( ...
                Na_k_2 + K_k_2 + Cl_k_2 + HCO3_k_2 - ...
                Na_s_2 - Cl_s_2 - K_s_2 - HCO3_s_2 + p.X_k ./ R_k_2);
            
            % Nernst potentials
            E_K_k_2 = p.R_g * p.T / (p.z_K * p.F) * log(K_s_2 ./ K_k_2);
            E_Na_k_2 = p.R_g * p.T / (p.z_Na * p.F) * log(Na_s_2 ./ Na_k_2);
            E_Cl_k_2 = p.R_g * p.T / (p.z_Cl * p.F) * log(Cl_s_2 ./ Cl_k_2);
            E_NBC_k_2 = p.R_g * p.T / (p.z_NBC * p.F) * ...
                log((Na_s_2 .* HCO3_s_2.^2) ./ (Na_k_2 .* HCO3_k_2.^2));
            E_BK_k_2 = p.R_g * p.T / (p.z_K * p.F) * log(K_p_2 ./ K_k_2);
            
            % Flux through the Sodium Potassium pump
            J_NaK_k_2 = p.J_NaK_max * Na_k_2.^1.5 ./ ...
                (Na_k_2.^1.5 + p.K_Na_k^1.5) .* ...
                K_s_2 ./ (K_s_2 + p.K_K_s);
            
            % Membrane voltage
            v_k_2 = (p.g_Na_k * E_Na_k_2 + p.g_K_k * E_K_k_2 + ...
                p.g_Cl_k * E_Cl_k_2 + p.g_NBC_k * E_NBC_k_2 + ...
                p.g_BK_k * w_k_2 .* E_BK_k_2 - ...
                J_NaK_k_2 * p.F / p.C_correction) ./ ...
                (p.g_Na_k + p.g_K_k + p.g_Cl_k + p.g_NBC_k + ...
                p.g_BK_k * w_k_2);
            
            % Fluxes
            J_BK_k_2 = p.g_BK_k / p.F * w_k_2 .* ...
                (v_k_2 - E_BK_k_2) * p.C_correction;
            J_K_k_2 = p.g_K_k / p.F * (v_k_2 - E_K_k_2) * p.C_correction;
            
            J_Na_k_2 = p.g_Na_k / p.F * (v_k_2 - E_Na_k_2) * p.C_correction;
            J_NBC_k_2 = p.g_NBC_k / p.F * (v_k_2 - E_NBC_k_2) * p.C_correction;
            J_KCC1_k_2 = self.flux_ft_2(t) .* p.g_KCC1_k / p.F * p.R_g * ...
                p.T / p.F .* ...
                log((K_s_2 .* Cl_s_2) ./ (K_k_2 .* Cl_k_2)) * p.C_correction;
            J_NKCC1_k_2 = self.flux_ft_2(t) * p.g_NKCC1_k / p.F * p.R_g * ...
                p.T / p.F .* log((Na_s_2 .* K_s_2 .* Cl_s_2.^2) ./ ...
                (Na_k_2 .* K_k_2 .* Cl_k_2.^2)) * p.C_correction;
            
            
            % Aditional equations
            w_inf_2 = 0.5 * ...
                (1 + tanh((v_k_2 + p.v_6) / p.v_4));
            
            phi_w_2 = p.psi_w * cosh((v_k_2 + p.v_6) / (2*p.v_4));
            
            %% Conservation Equations
            % Differential Equations in the Astrocyte
            du(idx.N_K_k_2, :) = -J_K_k_2 + 2*J_NaK_k_2 + J_NKCC1_k_2 + ...
                J_KCC1_k_2 - J_BK_k_2;
            du(idx.N_Na_k_2, :) = -J_Na_k_2 - 3*J_NaK_k_2 + J_NKCC1_k_2 + J_NBC_k_2;
            
            du(idx.N_HCO3_k_2, :) = 2*J_NBC_k_2;
            du(idx.N_Cl_k_2, :) = du(idx.N_Na_k_2, :) + du(idx.N_K_k_2, :) - ...
                du(idx.N_HCO3_k_2, :);
           
            du(idx.w_k_2, :) = phi_w_2 .* (w_inf_2 - w_k_2);
            
            % Differential Equations in the Perivascular space
            du(idx.K_p_2, :) = J_BK_k_2 ./ (R_k_2 * p.VR_pa) + J_KIR_i_2 ./ ...
                p.VR_ps - p.R_decay*(K_p_2 - p.K_p_min);
            
            % Differential Equations in the Synaptic Cleft
            du(idx.N_K_s_2, :) = p.k_C * self.input_f_2(t) - ...
                du(idx.N_K_k_2, :) - J_BK_k_2;
            du(idx.N_Na_s_2, :) = -p.k_C * self.input_f_2(t) - ...
                du(idx.N_Na_k_2, :);
            du(idx.N_HCO3_s_2, :) = -du(idx.N_HCO3_k_2, :);
            du = bsxfun(@times, self.enabled, du);
            if nargout == 2
               Uout = zeros(self.n_out, size(u, 2));
               Uout(self.idx_out.ft_2, :) = self.input_f_2(t);
               Uout(self.idx_out.v_k_2, :) = v_k_2;
               Uout(self.idx_out.K_s_2, :) = K_s_2;
               Uout(self.idx_out.K_p_2, :) = K_p_2;
               Uout(self.idx_out.J_BK_k_2, :) = J_BK_k_2;
               Uout(self.idx_out.w_inf_2, :) = w_inf_2;
               Uout(self.idx_out.phi_w_2, :) = phi_w_2;
              
               varargout = {Uout};
            end
        end        
        function K_p_2 = shared(self, ~, u)
            K_p_2 = u(self.index.K_p_2, :);
        end
        function f_2 = input_f_2(self, t)
            % The neuronal K+ input signal
            p = self.params;
            f_2 = zeros(size(t));
            
            ii = p.t_0 <= t & t < p.t_1;
            f_2(ii) = ...
                p.F_input * p.gab / ...
                (p.ga * p.gb) * ...
                (1 - (t(ii) - p.t_0) / p.delta_t).^(p.beta - 1) .* ...
                ((t(ii) - p.t_0) / p.delta_t).^(p.alpha - 1);
            f_2(p.t_2 <= t & t <= p.t_3) = -p.F_input;
        end
        function out_2 = flux_ft_2(self, t)
            % C_input Block function to switch channel on and off
            p = self.params;
            out_2 = ( ...
                0.5 * tanh((t - p.t_0) / 0.0005) - ...
                0.5 * tanh((t - p.t_1 - p.lengthpulse) / 0.0005));
            out_2 = out_2(:).';
        end
        function names = varnames(self)
            names = [fieldnames(self.index); fieldnames(self.idx_out)];
        end
    end    
end

function idx = indices()
% Index of parameters needing inital conditions 
idx.R_k_2 = 1;
idx.K_p_2 = 2;
idx.N_Na_k_2 = 3;
idx.N_K_k_2 = 4;
idx.N_Cl_k_2 = 5;
idx.N_HCO3_k_2 = 6;
idx.N_Na_s_2 = 7;
idx.N_K_s_2 = 8;
idx.N_HCO3_s_2 = 9;
idx.w_k_2 = 10;
end
function [idx, n] = output_indices()
% Index of all other output parameters
idx.ft_2 = 1;
idx.v_k_2 = 2;
idx.J_BK_k_2 = 3;
idx.K_s_2 = 4;
idx.K_p_2 = 5;
idx.w_inf_2 = 6;
idx.phi_w_2 = 7;
                    
n = numel(fieldnames(idx));
end
function params = parse_inputs(varargin)
parser = inputParser();
% Scaling Constants
parser.addParameter('L_p', 2.1e-9); % m uM^-1 s^-1
parser.addParameter('X_k', 12.41e-3); % uM m
parser.addParameter('R_tot', 8.79e-8); % m

% Input signal
parser.addParameter('startpulse', 200); % s     %%%%%%    
parser.addParameter('lengthpulse', 200); % s    %%%%%%
parser.addParameter('lengtht1', 10); % s
parser.addParameter('F_input', 2.5); % s
parser.addParameter('alpha', 2);% [-]
parser.addParameter('beta', 5);% [-]
parser.addParameter('delta_t', 10); % s

% Synpatic cleft
parser.addParameter('k_C', 7.35e-5); %uM m s^-1

% Perivascular space
parser.addParameter('VR_pa', 0.001);% [-]
parser.addParameter('VR_ps', 0.001);% [-]
parser.addParameter('R_decay', 0.05);% s^-1
parser.addParameter('K_p_min', 3e3);% uM

% Fluxes Constants
parser.addParameter('F', 9.65e4); %C mol^-1
parser.addParameter('R_g', 8.315); %J mol^-1 K^-1
parser.addParameter('T', 300); % K
parser.addParameter('g_K_k', 40); %mho m^-2
parser.addParameter('g_Na_k', 1.314); % mho m^-2
parser.addParameter('g_NBC_k', 7.57e-1); % mho m^-2
parser.addParameter('g_KCC1_k', 1e-2); % mho m^-2
parser.addParameter('g_NKCC1_k', 5.54e-2); % mho m^-2
parser.addParameter('J_NaK_max', 1.42e-3); % uM m s^-1
parser.addParameter('K_Na_k', 10000); % uM
parser.addParameter('K_K_s', 1500); % uM
parser.addParameter('G_BK_k', 4.3e3); % pS (later converted to mho m^-2)
parser.addParameter('A_ef_k', 3.7e-9); % m2
parser.addParameter('C_correction', 1e3); % [-]
parser.addParameter('J_max', 2880); %uM s^-1
parser.addParameter('K_act', 0.17); %uM
parser.addParameter('P_L', 0.0804); %uM
parser.addParameter('V_max', 20); %uM s^-1
parser.addParameter('k_pump', 0.24); %uM

% Additional Equations; Astrocyte Constants
parser.addParameter('g_Cl_k', 8.797e-1); % mho m^-2
parser.addParameter('z_K', 1);% [-]
parser.addParameter('z_Na', 1);% [-]
parser.addParameter('z_Cl', -1);% [-]
parser.addParameter('z_NBC', -1);% [-]
parser.addParameter('BK_end', 40);% [-]
parser.addParameter('K_ex', 0.26); %uM
parser.addParameter('B_ex', 11.35); %uM
parser.addParameter('K_G', 8.82); %uM
parser.addParameter('v_4', 14.5e-3); %V
parser.addParameter('v_5', 8e-3); %V
parser.addParameter('v_6', 22e-3); %V
parser.addParameter('psi_w', 2.664); %s^-1

parser.parse(varargin{:})
params = parser.Results;
params.g_BK_k = params.G_BK_k*1e-12 / params.A_ef_k;
params.t_0 = params.startpulse;
params.t_1 = params.t_0 + params.lengtht1;
params.t_2 = params.t_0 + params.lengthpulse;
params.t_3 = params.t_1 + params.lengthpulse;
params.gab = factorial(params.alpha + params.beta - 1);
params.ga = factorial(params.alpha - 1);
params.gb = factorial(params.beta - 1);
end
function u0 = initial_conditions(idx)
% Inital estimations of parameters from experimental data
u0 = zeros(length(fieldnames(idx)), 1);

u0(idx.R_k_2) = 0.061e-6;
u0(idx.N_Na_k_2) = 0.99796e-3;
u0(idx.N_K_k_2) = 5.52782e-3;
u0(idx.N_HCO3_k_2) = 0.58804e-3;
u0(idx.N_Cl_k_2) = 0.32879e-3;
u0(idx.N_Na_s_2) = 4.301041e-3;
u0(idx.N_K_s_2) = 0.0807e-3;
u0(idx.N_HCO3_s_2) = 0.432552e-3;
u0(idx.K_p_2) = 3e3;
u0(idx.w_k_2) = 0.1815e-3;
end
