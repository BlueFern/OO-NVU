classdef Astrocyte < handle
    % The 'Astrocyte' code contains the following sections of the model:
    %   The Neuron, Synaptic cleft, the Astrocyte and the Perivascular Space
    %   Currently there is no content under the Neuron sub-section
    %   Please refer to the relevient sections in the documentation for 
    %   full information on the equations and variable names.

    properties
        params
        u0
        index
        n_out
        idx_out
        enabled
    end
    methods
        function self = Astrocyte(varargin)
            self.params = parse_inputs(varargin{:});
            self.index = indices();
            self.u0 = initial_conditions(self.index);
            self.enabled = true(size(self.u0));
            [self.idx_out, self.n_out] = output_indices();
        end
        function [du, varargout] = rhs(self, t, u, J_KIR_i, R)
            
            % Initalise inputs and parameters
            t = t(:).';
            p = self.params;
            idx = self.index;
            
            R_k = u(idx.R_k, :);
            K_p = u(idx.K_p, :);
            Ca_p = u(idx.Ca_p, :);
            
            N_Na_k = u(idx.N_Na_k, :);
            N_K_k = u(idx.N_K_k, :);
            N_Cl_k = u(idx.N_Cl_k, :);
            N_HCO3_k = u(idx.N_HCO3_k, :);
            
            N_Na_s = u(idx.N_Na_s, :);
            N_K_s = u(idx.N_K_s, :);
            N_HCO3_s = u(idx.N_HCO3_s, :);
                   
            w_k = u(idx.w_k, :);
            
            i_k = u(idx.i_k, :);
            c_k = u(idx.c_k, :);
            h_k = u(idx.h_k, :);
            s_k = u(idx.s_k, :);
            z_k = u(idx.z_k, :); %TRPV
            eet_k = u(idx.eet_k, :);

            du = zeros(size(u));
           
            %% Synaptic Cleft: 
            % Electroneutrality condition 
            N_Cl_s = N_Na_s + N_K_s - N_HCO3_s;
            % Volume-surface ratio
            R_s = p.R_tot - R_k; 
            
            % Scale concentrations to get actual concentrations,
            K_s = N_K_s ./ R_s;
            Na_s = N_Na_s ./ R_s;
            Cl_s = N_Cl_s ./ R_s;
            HCO3_s = N_HCO3_s ./ R_s;
            Na_k = N_Na_k ./ R_k;
            K_k = N_K_k./R_k;
            Cl_k = N_Cl_k ./ R_k;
            HCO3_k = N_HCO3_k ./ R_k;
            
            %% Astrocyte 
            % Volume-surface Ratio; Scaling ODE
            du(idx.R_k, :) = p.L_p * ( ...
                Na_k + K_k + Cl_k + HCO3_k - ...
                Na_s - Cl_s - K_s - HCO3_s + p.X_k ./ R_k);
                                    
            % Nernst potentials
            E_K_k = p.R_g * p.T / (p.z_K * p.F) * log(K_s ./ K_k);
            E_Na_k = p.R_g * p.T / (p.z_Na * p.F) * log(Na_s ./ Na_k);
            E_Cl_k = p.R_g * p.T / (p.z_Cl * p.F) * log(Cl_s ./ Cl_k);
            E_NBC_k = p.R_g * p.T / (p.z_NBC * p.F) * ...
                log((Na_s .* HCO3_s.^2) ./ (Na_k .* HCO3_k.^2));
            E_BK_k = p.R_g * p.T / (p.z_K * p.F) * log(K_p ./ K_k);
            E_TRPV_k = p.R_g * p.T / (p.z_Ca * p.F) * log(Ca_p./c_k); %TRPV
            
            % Flux through the Sodium Potassium pump
            J_NaK_k = p.J_NaK_max * Na_k.^1.5 ./ ...
                (Na_k.^1.5 + p.K_Na_k^1.5) .* ...
                K_s ./ (K_s + p.K_K_s);
            
            % Membrane voltage
            v_k = (p.g_Na_k * E_Na_k + p.g_K_k * E_K_k + p.g_TRPV_k * z_k  .* E_TRPV_k + ...
               p.g_Cl_k * E_Cl_k +  p.g_NBC_k * E_NBC_k + ...
                p.g_BK_k * w_k .* E_BK_k - ...                 
               J_NaK_k * p.F / p.C_correction) ./ ...
                (p.g_Na_k +p.g_K_k +   p.g_Cl_k + p.g_NBC_k + p.g_TRPV_k * z_k + ... %TRPV
                p.g_BK_k * w_k);
            % 
            % Fluxes
            J_BK_k = p.g_BK_k / p.F * w_k .* ...
                (v_k - E_BK_k) * p.C_correction;
            J_K_k = p.g_K_k / p.F * (v_k - E_K_k) * p.C_correction;
            
            J_Na_k = p.g_Na_k / p.F * (v_k - E_Na_k) * p.C_correction;
            J_NBC_k = p.g_NBC_k / p.F * (v_k - E_NBC_k) * p.C_correction;
            J_KCC1_k = self.flux_ft(t) .* p.g_KCC1_k / p.F * p.R_g * ...
                p.T / p.F .* ...
                log((K_s .* Cl_s) ./ (K_k .* Cl_k)) * p.C_correction;
            J_NKCC1_k = self.flux_ft(t) * p.g_NKCC1_k / p.F * p.R_g * ...
                p.T / p.F .* log((Na_s .* K_s .* Cl_s.^2) ./ ...
                (Na_k .* K_k .* Cl_k.^2)) * p.C_correction;
             
            
            %% Calcium Equations
            % Flux
            J_IP3 = p.J_max * (...
                i_k ./ (i_k + p.K_I) .* ...
                c_k ./ (c_k + p.K_act) .* h_k).^3 .* (1 - c_k ./ s_k);
            J_ER_leak = p.P_L * (1 - c_k ./ s_k);
            J_pump = p.V_max * c_k.^2 ./ (c_k.^2 + p.k_pump^2);
            J_TRPV_k = -0.5* p.g_TRPV_k/p.F * z_k .*(v_k - E_TRPV_k) *... %TRPV4 E_TRPV_k
               p.C_correction ;
            
            % Other equations
            B_cyt = 1 ./ (1 + p.BK_end + p.K_ex * p.B_ex ./ ...
                (p.K_ex + c_k).^2);
            G = (self.input_rho(t) + p.delta) ./ ...
                (p.K_G + self.input_rho(t) + p.delta);
            v_3 = -p.v_5 / 2 * tanh((c_k - p.Ca_3) / p.Ca_4)+p.v_7;
            
            %% Parent Calcium equations
            w_inf = 0.5 * ...
                (1 + tanh((v_k + p.eet_shift * eet_k - v_3) / p.v_4));
            phi_w = p.psi_w * cosh((v_k - v_3) / (2*p.v_4));
            
             %% TRPV Channel open probabilty equations
            H_Ca_k = c_k./p.gam_cai_k+Ca_p./p.gam_cae_k;
            eta=(R-p.R_0_passive_k)./(p.R_0_passive_k);
            zinf_k=(1./(1+exp(-(eta-p.epshalf_k)./p.kappa_k))).*((1./(1+H_Ca_k)).*(H_Ca_k+tanh(((v_k)-p.v1_TRPV_k)./p.v2_TRPV_k))); %Define epsilon
            t_Ca_k= p.t_TRPV_k./Ca_p;
            
            %% Conservation Equations
            % Differential Equations in the Astrocyte
            du(idx.N_K_k, :) = (-J_K_k + 2*J_NaK_k + J_NKCC1_k + ...
                J_KCC1_k - J_BK_k);
            du(idx.N_Na_k, :) = -J_Na_k - 3*J_NaK_k + J_NKCC1_k + J_NBC_k;
            
            du(idx.N_HCO3_k, :) = 2*J_NBC_k;
            du(idx.N_Cl_k, :) = du(idx.N_Na_k, :) + du(idx.N_K_k, :) - ...
                du(idx.N_HCO3_k, :);
            % Differential Calcium Equations in Astrocyte
            du(idx.c_k, :) = B_cyt .* (J_IP3 - J_pump + J_ER_leak)+J_TRPV_k./R_k; %TRPV added
            du(idx.s_k, :) = -du(idx.c_k, :) ./(p.VR_ER_cyt) + (J_TRPV_k./(R_k.*p.VR_ER_cyt));
            du(idx.h_k, :) = p.k_on * (p.K_inh - (c_k + p.K_inh) .* h_k);
            du(idx.i_k, :) = p.r_h * G - p.k_deg * i_k;
            du(idx.z_k, :) = p.trpv_switch.*((zinf_k-z_k)./(t_Ca_k.*Ca_p)) ; %TRPV
            du(idx.eet_k, :) = p.V_eet * max(c_k - p.c_k_min, 0) - ...
                p.k_eet * eet_k;
            du(idx.w_k, :) =  phi_w .* (w_inf - w_k);
            
            % Differential Equations in the Perivascular space
            du(idx.K_p, :) = J_BK_k ./ (R_k*p.VR_pa) + J_KIR_i ./ ...
                p.VR_ps - p.R_decay*(K_p - p.K_p_min);
            du(idx.Ca_p, :) = 0; %- J_TRPV_k ./ (p.VR_pa) - p.Ca_decay_k*(Ca_p - p.Capmin_k);
            % Differential Equations in the Synaptic Cleft
            du(idx.N_K_s, :) = p.k_C * self.input_f(t) - ...
                du(idx.N_K_k, :) - J_BK_k;
            du(idx.N_Na_s, :) = -p.k_C * self.input_f(t) - ...
                du(idx.N_Na_k, :);
            du(idx.N_HCO3_s, :) = -du(idx.N_HCO3_k, :);
            du = bsxfun(@times, self.enabled, du);
            if nargout == 2
               Uout = zeros(self.n_out, size(u, 2));
               Uout(self.idx_out.ft, :) = self.input_f(t);
               Uout(self.idx_out.v_k, :) = v_k;
               Uout(self.idx_out.K_s, :) = K_s;
               Uout(self.idx_out.K_p, :) = K_p;
               Uout(self.idx_out.J_BK_k, :) = J_BK_k;
               Uout(self.idx_out.rho, :) = self.input_rho(t);
               Uout(self.idx_out.B_cyt, :) = B_cyt;               
               Uout(self.idx_out.G, :) = G;               
               Uout(self.idx_out.v_3, :) = v_3;
               Uout(self.idx_out.w_inf, :) = w_inf;
               Uout(self.idx_out.phi_w, :) = phi_w;
               Uout(self.idx_out.J_IP3, :) = J_IP3;
               Uout(self.idx_out.J_pump, :) = J_pump;
               Uout(self.idx_out.J_ER_leak, :) = J_ER_leak;
               Uout(self.idx_out.J_TRPV_k, :) = J_TRPV_k;
%                Uout(self.idx_out.z_k, :) = z_k;
%                Uout(self.idx_out.w_k, :) = w_k;
               Uout(self.idx_out.E_BK_k, :) =  E_BK_k;
               Uout(self.idx_out.E_TRPV_k, :) =  E_TRPV_k;
               varargout = {Uout};
            end
        end        
        function K_p = shared(self, ~, u)
            K_p = u(self.index.K_p, :);
        end
        function f = input_f(self, t)
            % The neuronal K+ input signal
            p = self.params;
            f = zeros(size(t));
            ii = p.t_0 <= t & t < p.t_1;
            f(ii) = p.input_switch .*(...
                p.F_input * p.gab / ...
                (p.ga * p.gb) * ...
                (1 - (t(ii) - p.t_0) / p.delta_t).^(p.beta - 1) .* ...
                ((t(ii) - p.t_0) / p.delta_t).^(p.alpha - 1));
            f(p.t_2 <= t & t <= p.t_3) = p.input_switch .*(-p.F_input);
        end
        function rho = input_rho(self, t)
            % Input signal; the smooth pulse function rho
            p = self.params;
            rho = p.input_switch.*((p.Amp - p.base) * ( ...
                0.5 * tanh((t - p.t_0) / p.theta_L) - ...
               0.5 * tanh((t - p.t_2) / p.theta_R))) + p.base;
        end
        function out = flux_ft(self, t)
            % C_input Block function to switch channel on and off
            p = self.params;
            out = p.input_switch .*(( ...
                0.5 * tanh((t - p.t_0) / 0.0005) - ...
                0.5 * tanh((t - p.t_1 - p.lengthpulse) / 0.0005)));
            out = out(:).';
        end
        function names = varnames(self)
            names = [fieldnames(self.index); fieldnames(self.idx_out)];
        end
    end    
end

function idx = indices()
% Index of parameters needing inital conditions 
idx.R_k = 1;
idx.K_p = 2;
idx.N_Na_k = 3;
idx.N_K_k = 4;
idx.N_Cl_k = 5;
idx.N_HCO3_k = 6;
idx.N_Na_s = 7;
idx.N_K_s = 8;
idx.N_HCO3_s = 9;
idx.w_k = 10;
idx.i_k = 11;
idx.c_k = 12;
idx.h_k = 13;
idx.s_k = 14;
idx.eet_k = 15;
idx.z_k = 16;
idx.Ca_p = 17;
end
function [idx, n] = output_indices()
% Index of all other output parameters
idx.ft = 1;
idx.v_k = 2;
idx.J_BK_k = 3;
idx.K_s = 4;
idx.K_p = 5;
idx.rho = 6;
idx.B_cyt = 7;
idx.G = 8;
idx.v_3 = 9;
idx.w_inf = 10;
idx.phi_w = 11;
idx.J_IP3 = 12;
idx.J_pump = 13;
idx.J_ER_leak = 14;
idx.J_TRPV_k = 15;
% idx.z_k = 16;
% idx.w_k = 17;
idx.E_BK_k = 16;
idx.E_TRPV_k = 17;
n = numel(fieldnames(idx));
end
function params = parse_inputs(varargin)
parser = inputParser();
% Scaling Constants
parser.addParameter('L_p', 2.1e-9); % m uM^-1 s^-1
parser.addParameter('X_k', 12.41e-3); % uM m
parser.addParameter('R_tot', 8.79e-8); % m

% Input signal
parser.addParameter('startpulse', 200); % s
parser.addParameter('lengthpulse', 200); % s
parser.addParameter('lengtht1', 10); % s
parser.addParameter('F_input', 2.5); % s
parser.addParameter('alpha', 2);% [-]
parser.addParameter('beta', 5);% [-]
parser.addParameter('delta_t', 10); % s
parser.addParameter('input_switch', 1); % 

% Calcium in the Astrocyte Equations Constants
parser.addParameter('Amp', 0.7);
parser.addParameter('base', 0.1);
parser.addParameter('theta_L', 1);
parser.addParameter('theta_R', 1);
parser.addParameter('delta', 1.235e-2); %THIS MAY BE WRONG

parser.addParameter('VR_ER_cyt', 0.185)
parser.addParameter('k_on', 2); %uM s^-1
parser.addParameter('K_inh', 0.1); %uM
parser.addParameter('r_h', 4.8); % uM
parser.addParameter('k_deg', 1.25); % s^-1
parser.addParameter('V_eet', 72); % uM
parser.addParameter('k_eet', 7.2); % uM
parser.addParameter('c_k_min', 0.1); % uM

parser.addParameter('v_7', -15e-3); %V

parser.addParameter('Ca_3', 0.4);
parser.addParameter('Ca_4', 0.15);
parser.addParameter('eet_shift', 2e-3);

parser.addParameter('K_I', 0.03); % uM
%TRPV4
parser.addParameter('Capmin_k', 5); %uM

parser.addParameter('C_astr_k', 40);%pF
parser.addParameter('gamma_k', 834.3);%mV/uM
parser.addParameter('gam_cae_k', 200); %uM
parser.addParameter('gam_cai_k', 0.01); %uM
parser.addParameter('epshalf_k', 0.1);
parser.addParameter('kappa_k', 0.1);
parser.addParameter('v1_TRPV_k', 0.120); %mV
parser.addParameter('v2_TRPV_k', 0.013); %mV
parser.addParameter('t_TRPV_k', 0.9); %mV
parser.addParameter('R_0_passive_k', 20e-6)
parser.addParameter('trpv_switch', 1)
parser.addParameter('Ca_decay_k', 0.5)

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
parser.addParameter('z_Ca', 2);% [-]

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
params.g_TRPV_k = (50 * 1e-12)/params.A_ef_k;%mho m^-2
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

u0(idx.R_k) = 0.061e-6;
u0(idx.N_Na_k) = 0.99796e-3;
u0(idx.N_K_k) = 5.52782e-3;
u0(idx.N_HCO3_k) = 0.58804e-3;
u0(idx.N_Cl_k) = 0.32879e-3;
u0(idx.N_Na_s) = 4.301041e-3;
u0(idx.N_K_s) = 0.0807e-3;
u0(idx.N_HCO3_s) = 0.432552e-3;
u0(idx.K_p) = 3e3;
u0(idx.w_k) = 0.1815e-3;
u0(idx.c_k) = 0.05e-3;
u0(idx.s_k) = 400;
u0(idx.h_k) = 0.1e-3;
u0(idx.i_k) = 0.01e-3;
u0(idx.eet_k) = 0.1e-3;
u0(idx.z_k) = 0;
u0(idx.Ca_p) = 5; %5.1
end
