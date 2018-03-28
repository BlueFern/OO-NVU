classdef Astrocyte < handle
    % The 'Astrocyte' code contains the following sections of the model:
    % The Neuron, Synaptic cleft, the Astrocyte and the Perivascular Space
    % Currently there is no content under the Neuron sub-section
    % Please refer to the relevient sections in the documentation for
    % full information on the equations and variable names.
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
            self.index = indices(self);
            self.u0 = initial_conditions(self.index,self);
            self.enabled = true(size(self.u0));
            [self.idx_out, self.n_out] = output_indices(self);
        end

        function [du, varargout] = rhs(self, t, u, J_KIR_i, R, J_VOCC_i, NO_n, NO_i, J_K_NEtoSC, J_Na_NEtoSC, Glu, ATPg)
            % Initalise inputs and parameters
            t = t(:).';
            p = self.params;
            idx = self.index;
            
            v_k = u(idx.v_k, :);
            
            K_p = u(idx.K_p, :);
            Ca_p = u(idx.Ca_p, :);
            Na_k = u(idx.Na_k, :);
            K_k = u(idx.K_k, :);
            Cl_k = u(idx.Cl_k, :);
            HCO3_k = u(idx.HCO3_k, :);
            Na_s = u(idx.Na_s, :);
            K_s = u(idx.K_s, :);
            HCO3_s = u(idx.HCO3_s, :);
            w_k = u(idx.w_k, :);
            I_k = u(idx.I_k, :);
            
            Ca_k = u(idx.Ca_k, :);
            h_k = u(idx.h_k, :);
            s_k = u(idx.s_k, :);
            m_k = u(idx.m_k, :); 
            eet_k = u(idx.eet_k, :);
            NO_k = u(idx.NO_k, :);
            du = zeros(size(u));
            
            %% Synaptic Cleft:

            % Electroneutrality condition
            Cl_s = Na_s + K_s - HCO3_s;

            % Volume ratio of SC to AC
            VR_sa = p.R_s/p.R_k;
            
            % Input of K+ and Na+ to the SC (assuming that the SC is a small part of the ECS and everything that happens to the ECS also happens to the SC)
            J_K_NEtoSC_k = J_K_NEtoSC * 1e3; % Convert from mM/s to uM/s
            J_Na_NEtoSC_k = J_Na_NEtoSC * 1e3; % Convert from mM/s to uM/s
            
            %% Astrocyte

            % Nernst potentials (in mV)
            E_K_k = p.ph / p.z_K * log(K_s ./ K_k);
            E_Na_k = p.ph / p.z_Na * log(Na_s ./ Na_k);
            E_Cl_k = p.ph / p.z_Cl * log(Cl_s ./ Cl_k);
            E_NBC_k = p.ph / p.z_NBC * log((Na_s .* HCO3_s.^2) ./ (Na_k .* HCO3_k.^2));
            E_BK_k = p.ph / p.z_K * log(K_p ./ K_k);        
            E_TRPV_k = p.ph / p.z_Ca * log(Ca_p ./ Ca_k);                                     % Nernst potential TRPV
            
            % Flux through the Sodium Potassium pump
            J_NaK_k = p.J_NaK_max * Na_k.^1.5 ./ (Na_k.^1.5 + p.K_Na_k^1.5) .* K_s ./ (K_s + p.K_K_s);
            
            % Fluxes
            J_BK_k = p.G_BK_k * w_k .* (v_k - E_BK_k);      % BK flux (uM/s)
            J_BK_p = J_BK_k ./ p.VR_pa;                         % K+ influx into the PVS (uM/s)
            J_K_k = p.G_K_k * (v_k - E_K_k);
            J_Na_k = p.G_Na_k * (v_k - E_Na_k);
            J_NBC_k = p.G_NBC_k * (v_k - E_NBC_k);
            J_KCC1_k = p.G_KCC1_k * p.ph .* log((K_s .* Cl_s) ./ (K_k .* Cl_k));
            J_NKCC1_k = p.G_NKCC1_k * p.ph .* log((Na_s .* K_s .* Cl_s.^2) ./ (Na_k .* K_k .* Cl_k.^2));
            J_Cl_k = p.G_Cl_k * (v_k - E_Cl_k);
            
            %% Calcium Equations
            % Flux
            J_IP3 = p.J_max * ( I_k ./ (I_k + p.K_I) .*  Ca_k ./ (Ca_k + p.K_act) .* h_k).^3 .* (1 - Ca_k ./ s_k);
            J_ER_leak = p.P_L * (1 - Ca_k ./ s_k);
            J_pump = p.V_max * Ca_k.^2 ./ (Ca_k.^2 + p.k_pump^2);
%             I_TRPV_k = p.G_TRPV_k * m_k .* (v_k - E_TRPV_k); % TRPV4 current
%             J_TRPV_k = - I_TRPV_k / (p.z_Ca * p.C_astr_k * p.gamma_i); % TRPV4 flux [uM/s]
            
            J_TRPV_k = p.G_TRPV_k * m_k .* (v_k - E_TRPV_k);

            rho = p.rho_min + (p.rho_max - p.rho_min)/p.Glu_max * Glu;
            
            % Other equations
            B_cyt = 1 ./ (1 + p.BK_end + p.K_ex * p.B_ex ./ (p.K_ex + Ca_k).^2);
            G = (rho + p.delta) ./ (p.K_G + rho + p.delta);
            v_3 = p.v_6 - p.v_5 / 2 * tanh((Ca_k - p.Ca_3) / p.Ca_4);
            
            % Parent Calcium equations
            
            w_inf = 0.5 * (1 + tanh((v_k + p.eet_shift * eet_k - v_3) / p.v_4));
            phi_w = p.psi_w * cosh((v_k - v_3) / (2 * p.v_4));
		   
            %% TRPV Channel open probability equations
            H_Ca_k = Ca_k ./ p.gam_cai_k + Ca_p ./ p.gam_cae_k;
            eta = (R - p.R_init) ./ (p.R_init);
            minf_k = (1 ./ (1 + exp(-(eta - p.epshalf_k) ./ p.kappa_k))) .* ((1 ./ (1 + H_Ca_k)) .* (H_Ca_k + tanh((v_k - p.v1_TRPV_k) ./ p.v2_TRPV_k))); 
            
            %% NO pathway
            tau_nk = p.x_nk ^ 2 ./  (2 * p.D_cNO);
            tau_ki = p.x_ki ^ 2 ./  (2 * p.D_cNO);
            p_NO_k = 0;
            c_NO_k = p.k_O2_k * NO_k.^2 * p.O2_k; % [uM/s]
            d_NO_k = (NO_n - NO_k) ./ tau_nk + (NO_i - NO_k) ./ tau_ki;


            %% Conservation Equations
            
            du(idx.v_k, :)    = p.gamma_i .* ( -J_BK_k - J_K_k - J_Cl_k - J_NBC_k - J_Na_k - J_NaK_k - 2*J_TRPV_k);
            
            % Differential Equations in the Astrocyte
            du(idx.K_k, :)    = -J_K_k + 2*J_NaK_k + J_NKCC1_k + J_KCC1_k - J_BK_k;
            du(idx.Na_k, :)   = -J_Na_k - 3*J_NaK_k + J_NKCC1_k + J_NBC_k;
            du(idx.HCO3_k, :) = 2*J_NBC_k;
            du(idx.Cl_k, :)   = du(idx.Na_k, :) + du(idx.K_k, :) - du(idx.HCO3_k, :) + 2*du(idx.Ca_k, :);
            
            % Differential Calcium Equations in Astrocyte
            du(idx.Ca_k, :)     = B_cyt .* (J_IP3 - J_pump + J_ER_leak - J_TRPV_k/p.r_buff);           
            du(idx.s_k, :)      = -(B_cyt .* (J_IP3 - J_pump + J_ER_leak)) ./ (p.VR_ER_cyt);
            du(idx.h_k, :)      = p.k_on * (p.K_inh - (Ca_k + p.K_inh) .* h_k);
            du(idx.I_k, :)      = p.r_h * G - p.k_deg * I_k;
            du(idx.m_k, :)      = p.trpv_switch .* ((minf_k - m_k) ./ p.t_TRPV_k) ; % TRPV open probability, the trpv_switch can switch the turn the channel on or off.
            du(idx.eet_k, :)    = p.V_eet * max(Ca_k - p.Ca_k_min, 0) - p.k_eet * eet_k;
            du(idx.w_k, :)      = phi_w .* (w_inf - w_k);
            
            % Differential Equations in the Perivascular space
            du(idx.K_p, :)      = J_BK_k ./ (p.VR_pa) + J_KIR_i ./ p.VR_ps - p.R_decay * (K_p - p.K_p_min) ;
            du(idx.Ca_p, :)     = J_TRPV_k./ p.VR_pa + J_VOCC_i ./ p.VR_ps - p.Ca_decay_k .* (Ca_p - p.Capmin_k); % calcium concentration in PVS
           
            % Differential Equations in the Synaptic Cleft
            du(idx.K_s, :)    = 1./VR_sa * (J_K_k - 2 * J_NaK_k - J_NKCC1_k - J_KCC1_k) + J_K_NEtoSC_k;
            du(idx.Na_s, :)   = 1./VR_sa * (J_Na_k + 3*J_NaK_k - J_NKCC1_k - J_NBC_k) + J_Na_NEtoSC_k;
            du(idx.HCO3_s, :) = 1./VR_sa * (-2*J_NBC_k);
 
            % NO pathway
            du(idx.NO_k, :) = p_NO_k - c_NO_k + d_NO_k;
             
            du = bsxfun(@times, self.enabled, du);
            if nargout == 2

                Uout = zeros(self.n_out, size(u, 2));
                Uout(self.idx_out.J_K_NEtoSC, :) = J_K_NEtoSC;
                Uout(self.idx_out.K_s, :) = K_s;
                Uout(self.idx_out.K_p, :) = K_p;
                Uout(self.idx_out.J_BK_k, :) = J_BK_k;
                Uout(self.idx_out.rho, :) = rho;
                Uout(self.idx_out.B_cyt, :) = B_cyt;
                Uout(self.idx_out.G, :) = G;
                Uout(self.idx_out.v_3, :) = v_3;
                Uout(self.idx_out.J_IP3, :) = J_IP3;
                Uout(self.idx_out.J_pump, :) = J_pump;
                Uout(self.idx_out.J_ER_leak, :) = J_ER_leak;
                Uout(self.idx_out.J_TRPV_k, :) = J_TRPV_k;
                Uout(self.idx_out.E_BK_k, :) = E_BK_k;
                Uout(self.idx_out.J_NBC_k, :) = J_NBC_k;
                Uout(self.idx_out.E_Na_k, :) = E_Na_k;
                Uout(self.idx_out.E_K_k, :) = E_K_k;
                Uout(self.idx_out.E_NBC_k, :) = E_NBC_k;
                Uout(self.idx_out.E_Cl_k, :) = E_Cl_k;
%                 Uout(self.idx_out.I_TRPV_k, :) = I_TRPV_k;
                Uout(self.idx_out.J_K_k, :) = J_K_k;
                Uout(self.idx_out.J_Na_k, :) = J_Na_k;
                Uout(self.idx_out.K_k, :) = K_k;
                Uout(self.idx_out.w_inf, :) = w_inf;
                Uout(self.idx_out.phi_w, :) = phi_w;
                Uout(self.idx_out.J_BK_p, :) = J_BK_p;
                Uout(self.idx_out.E_TRPV_k, :) = E_TRPV_k;
                Uout(self.idx_out.Na_k, :) = Na_k;
                Uout(self.idx_out.J_KCC1_k, :) = J_KCC1_k;
                Uout(self.idx_out.J_NKCC1_k, :) = J_NKCC1_k;
                Uout(self.idx_out.J_NaK_k, :) = J_NaK_k;
                Uout(self.idx_out.J_Cl_k, :) = J_Cl_k;
                
                
                
                varargout = {Uout};
            end
        end
        function [K_p, NO_k, Na_k] = shared(self, ~, u)
            idx = self.index;
            K_p = u(self.index.K_p, :);
            NO_k = u(idx.NO_k, :);
            Na_k = u(idx.Na_k, :);
            
        end
        
        function names = varnames(self)
            names = [fieldnames(self.index); fieldnames(self.idx_out)];
        end
    end   
end

function idx = indices(self)
% Index of state variables
    idx.NO_k = 1;
    idx.K_p = 2;
    idx.Na_k = 3;
    idx.K_k = 4;
    idx.Cl_k = 5;
    idx.HCO3_k = 6;
    idx.Na_s = 7;
    idx.K_s = 8;
    idx.HCO3_s = 9;
    idx.w_k = 10;
    idx.I_k = 11;
    idx.Ca_k = 12;
    idx.h_k = 13;
    idx.s_k = 14;
    idx.eet_k = 15;
    idx.m_k = 16;
    idx.Ca_p = 17;
    idx.v_k = 18;
end

function [idx, n] = output_indices(self)
    % Index of all other output parameters
    idx.K_k  = 1;
    idx.v_k  = 2;
    idx.K_s  = 3;
    idx.K_p  = 4;
    idx.J_BK_k  = 5;
    idx.rho  = 6;
    idx.B_cyt  = 7;
    idx.G  = 8;
    idx.v_3  = 9;
    idx.w_inf  = 10;
    idx.phi_w  = 11;
    idx.J_IP3  = 12;
    idx.J_pump  = 13;
    idx.J_ER_leak  = 14;
    idx.J_TRPV_k  = 15;
    idx.E_BK_k  = 16;
    idx.J_NBC_k  = 17;
    idx.E_Na_k  = 18; 
    idx.E_K_k  = 19; 
    idx.E_NBC_k  = 20; 
    idx.E_Cl_k  = 21; 
%     idx.I_TRPV_k  = 22; 
    idx.J_NKCC1_k = 23;
    idx.J_K_k  = 24; 
    idx.J_Na_k  = 25;
    idx.J_BK_p = 26;
    idx.J_NaK_k = 27;
    idx.E_TRPV_k = 28;
    idx.Cak = 29;
    idx.Na_k = 30;
    idx.J_K_NEtoSC = 31;
    idx.J_KCC1_k = 32;
    idx.J_Cl_k = 33;
    
    n = numel(fieldnames(idx));
end

function params = parse_inputs(varargin)
    parser = inputParser();
    
    parser.addParameter('gamma_i', 1970); % [mV/uM] 
       
    % Conductances - now in uM mV^-1 s^-1 so consistent with other
    % fluxes, original conductances in mho m^-2 are commented
    % Converted by dividing by R_k=6e-8 and F=9.65e4
    parser.addParameter('G_BK_k', 10.25)        % g_BK_k = 225 * 1e-12 / 3.7e-9 = 6.08e-2;
    parser.addParameter('G_K_k', 6907.77);      % 40
    parser.addParameter('G_Na_k', 226.94);      % 1.314
    parser.addParameter('G_NBC_k', 130.74);     % 7.57e-1
    parser.addParameter('G_KCC1_k', 1.728);     % 1e-2
    parser.addParameter('G_NKCC1_k', 9.568);    % 5.54e-2
    parser.addParameter('G_Cl_k', 151.93);      % 8.797e-1

    parser.addParameter('J_NaK_max', 2.3667e4); % uM s^-1      1.42e-3
    
    % ECS K+ input parameters
    parser.addParameter('ECS_input', 9); 
    parser.addParameter('t0_ECS', 10000);
    parser.addParameter('tend_ECS', 20000);
    
    % Switches to turn on and off some things
    parser.addParameter('rhoSwitch', 1); 
    parser.addParameter('GluSwitch', 1); 
    parser.addParameter('trpv_switch', 1); 
    
    parser.addParameter('rho_min', 0.1);     
    parser.addParameter('rho_max', 0.7);  
    parser.addParameter('Glu_max', 1846);       % microM (one vesicle, Santucci2008)
    
    parser.addParameter('v_4', 8); %mV
    parser.addParameter('v_5', 15); %mV
    parser.addParameter('v_6', -55); %mV
    parser.addParameter('Ca_3', 0.4); % uM 
    parser.addParameter('Ca_4', 0.35); % uM 

    % Scaling Constants
    parser.addParameter('R_s', 2.79e-8); % m
    parser.addParameter('R_k', 6e-8); % m
    parser.addParameter('R_tot', 8.79e-8); % m
    
    % Calcium in the Astrocyte Equations Constants
    parser.addParameter('delta', 1.235e-2); 
    parser.addParameter('VR_ER_cyt', 0.185)
    parser.addParameter('k_on', 2); %uM s^-1
    parser.addParameter('K_inh', 0.1); %uM
    parser.addParameter('r_h', 4.8); % uM/s
    parser.addParameter('k_deg', 1.25); % s^-1
    parser.addParameter('V_eet', 72); % s^-1
    parser.addParameter('k_eet', 7.2); % s^-1
    parser.addParameter('Ca_k_min', 0.1); % uM
    parser.addParameter('eet_shift', 2); %mV/uM
    parser.addParameter('K_I', 0.03); % uM

    %TRPV4
    parser.addParameter('Capmin_k', 2000); %uM
    parser.addParameter('C_astr_k', 40);%pF
    parser.addParameter('gam_cae_k', 200); %uM
    parser.addParameter('gam_cai_k', 0.01); %uM
     parser.addParameter('epshalf_k', 0.1); % 
    parser.addParameter('kappa_k', 0.1);
    parser.addParameter('v1_TRPV_k', 120); %mV
    parser.addParameter('v2_TRPV_k', 13); %mV
    parser.addParameter('t_TRPV_k', 0.9); 
    parser.addParameter('R_init', 20); 
    parser.addParameter('Ca_decay_k', 0.5);
    parser.addParameter('G_TRPV_k', 3.15e-4); 
    parser.addParameter('r_buff', 0.05); % Rate at which Ca2+ from the TRPV4 channel at the endfoot is buffered compared to rest of channels on the astrocyte body [-]
    
    % Perivascular space
    parser.addParameter('VR_pa', 0.001);% [-]
    parser.addParameter('VR_ps', 0.001);% [-]
    parser.addParameter('R_decay', 0.15);% s^-1
    parser.addParameter('K_p_min', 3e3);% uM

    % Fluxes Constants
    parser.addParameter('F', 9.65e4); %C mol^-1; Faraday's constant
    parser.addParameter('R_g', 8.315); %J mol^-1 K^-1; Gas constant
    parser.addParameter('T', 300); % K; Temperature
        
    parser.addParameter('K_Na_k', 10000); % uM
    parser.addParameter('K_K_s', 1500); % uM
    
    parser.addParameter('J_max', 2880); %uM s^-1
    parser.addParameter('K_act', 0.17); %uM
    parser.addParameter('P_L', 0.0804); %uM s^-1
    parser.addParameter('V_max', 20); %uM s^-1
    parser.addParameter('k_pump', 0.24); %uM

    % Additional Equations; Astrocyte Constants
    parser.addParameter('z_K', 1);% [-]
    parser.addParameter('z_Na', 1);% [-]
    parser.addParameter('z_Cl', -1);% [-]
    parser.addParameter('z_NBC', -1);% [-]
    parser.addParameter('z_Ca', 2);% [-]
    parser.addParameter('BK_end', 40);% [-]
    parser.addParameter('K_ex', 0.26); %uM
    parser.addParameter('B_ex', 11.35); %uM
    parser.addParameter('K_G', 8.82); 
    parser.addParameter('psi_w', 2.664); %s^-1
    parser.addParameter('ph', 26.6995);         % RT/Farad where Farad is [C/mmol], also used in Neuron.m

    % NO Pathway
    parser.addParameter('D_cNO', 3300);         % [um^2 s^-1] ; Diffusion coefficient NO (Malinski1993)
    parser.addParameter('x_nk', 25);            % [um] ;  (M.E.)
    parser.addParameter('x_ki', 25);            % [um] ;  (M.E.)
    parser.addParameter('k_O2_k', 9.6e-6);      % [uM^-2 s^-1] ;  (Kavdia2002)
    parser.addParameter('O2_k', 200);           % [uM] ;  (M.E.)


    
    parser.parse(varargin{:})
    params = parser.Results;
    
end

function u0 = initial_conditions(idx,self)
    % Inital estimations of parameters from experimental data
    p = self.params;
    u0 = zeros(length(fieldnames(idx)), 1);
    u0(idx.Na_k) = 18268;   % [uM]
    u0(idx.K_k) = 92708;    % [uM]  
    u0(idx.HCO3_k) = 9131;  % [uM]
    u0(idx.Cl_k) = 7733;    % [uM]
    u0(idx.Na_s) = 150255;  % [uM]
    u0(idx.K_s) = 2837;     % [uM]
    u0(idx.HCO3_s) = 16881; % [uM]
    u0(idx.K_p) = 3045.1;   % [uM]
    u0(idx.w_k) = 1.703e-4; % [-]
    u0(idx.Ca_k) = 0.1612;  % [uM]
    u0(idx.s_k) = 480.8;    % [uM] 
    u0(idx.h_k) = 0.3828;   % [-]
    u0(idx.I_k) = 0.048299; % [uM]
    u0(idx.eet_k) = 0.6123; % [uM]
    u0(idx.m_k) = 0.5710;   % [-]
    u0(idx.Ca_p) = 1746.4;  % [uM]
    u0(idx.NO_k) = 0.1106;  % [uM]
    u0(idx.v_k) = -88.9;      % [mV] 
end
