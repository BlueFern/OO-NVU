classdef Neuron < handle
    %Contains neuron (soma/axon and dendrite) and extracellular space

    properties
        params
        u0
        index
        n_out
        idx_out
        enabled
    end
    methods
        function self = Neuron(varargin)
            self.params = parse_inputs(varargin{:});
            self.index = indices();
            self.u0 = initial_conditions(self.index);
            self.enabled = true(size(self.u0));
            [self.idx_out, self.n_out] = output_indices();
        end
        function [du, varargout] = rhs(self, t, u, NO_k, R, K_s, Na_s)
            t = t(:).';
            p = self.params;
            idx = self.index;
            
        %% State variables
            %% NO pathway 
            Ca_n = u(idx.Ca_n, :);                  
            nNOS_act_n = u(idx.nNOS_act_n, :);
            NO_n = u(idx.NO_n, :);
            
            %% Elshin neuron model
            v_sa = u(idx.v_sa, :);      % membrane potential of soma/axon, mV
            v_d = u(idx.v_d, :);        % membrane potential of dendrite, mV
            K_sa = u(idx.K_sa, :);      % K+ concentration of soma/axon, mM
            Na_sa = u(idx.Na_sa, :);    % Na+ concentration of soma/axon, mM
            Cl_sa = u(idx.Cl_sa, :);    % CL- concentration of soma/axon, mM
            K_d = u(idx.K_d, :);        % K+ concentration of dendrite, mM
            Na_d = u(idx.Na_d, :);      % Na+ concentration of dendrite, mM
            Cl_d = u(idx.Cl_d, :);      % Cl- concentration of dendrite, mM
            K_e = u(idx.K_e, :);        % K+ concentration of ECS, mM
            Na_e = u(idx.Na_e, :);      % Na+ concentration of ECS, mM
            Cl_e = u(idx.Cl_e, :);      % Cl- concentration of ECS, mM
            
            Buff = u(idx.Buff, :);      % Buffer concentration for K+ buffering in ECS, mM
            O2 = u(idx.O2, :);          % Tissue oxygen concentration, mM
            
            % Gating variables
            m1 = u(idx.m1, :);          % Activation gating variable, soma/axon NaP channel (Na+)
            m2 = u(idx.m2, :);          % Activation gating variable, soma/axon KDR channel (K+)
            m3 = u(idx.m3, :);          % Activation gating variable, soma/axon KA channel (K+)
            m4 = u(idx.m4, :);          % Activation gating variable, dendrite NaP channel (Na+)
            m5 = u(idx.m5, :);          % Activation gating variable, dendrite NMDA channel (Na+)
            m6 = u(idx.m6, :);          % Activation gating variable, dendrite KDR channel (K+)
            m7 = u(idx.m7, :);          % Activation gating variable, dendrite KA channel (K+)
            m8 = u(idx.m8, :);          % Activation gating variable, soma/axon NaT channel (Na+)
            
            h1 = u(idx.h1, :);          % Inactivation gating variable, soma/axon NaP channel (Na+)
            h2 = u(idx.h2, :);          % Inactivation gating variable, soma/axon KA channel (K+)
            h3 = u(idx.h3, :);          % Inactivation gating variable, dendrite NaP channel (Na+)
            h4 = u(idx.h4, :);          % Inactivation gating variable, dendrite NMDA channel (Na+)
            h5 = u(idx.h5, :);          % Inactivation gating variable, dendrite KA channel (K+)
            h6 = u(idx.h6, :);          % Inactivation gating variable, soma/axon NaT channel (Na+)
            
            
        %% Algebraic expressions       
            %% NO pathway
            w_NR2A      = self.input_Glu(t, v_sa) ./ (p.K_mA + self.input_Glu(t, v_sa)); %[-] 
            w_NR2B      = self.input_Glu(t, v_sa) ./ (p.K_mB + self.input_Glu(t, v_sa)); %[-]
            I_Ca        = (-4 * p.v_n * p.G_M * p.P_Ca_P_M * (p.Ca_ex / p.M)) ./ (1 + exp(-80 * (p.v_n + 0.02)))...         
                            .* (exp(2 * p.v_n * p.F / (p.R_gas * p.T))) ./ (1 - exp(2 * p.v_n * p.F / (p.R_gas * p.T))); %[fA]
%             I_Ca        = (-4 * v_sa/1000 * p.G_M * p.P_Ca_P_M * (p.Ca_ex / p.M)) ./ (1 + exp(-80 * (v_sa/1000 + 0.02)))...         
%                             .* (exp(2 * v_sa/1000 * p.F / (p.R_gas * p.T))) ./ (1 - exp(2 * v_sa/1000 * p.F / (p.R_gas * p.T))); %[fA]
            I_Ca_tot    = I_Ca .* (p.n_NR2A * w_NR2A + p.n_NR2B * w_NR2B); %[fA]
            phi_N       = 1 + p.Q1 * Ca_n + p.Q1 * p.Q2 * Ca_n.^2 + p.Q1 * p.Q2 * p.Q3 * Ca_n.^3 + p.Q1 * p.Q2 * p.Q3 * p.Q4 * Ca_n.^4; %[-]
            dphi_N      = p.Q1 + 2 * p.Q1 * p.Q2 * Ca_n + 3 * p.Q1 * p.Q2 * p.Q3 * Ca_n.^2 + 4 * p.Q1 * p.Q2 * p.Q3 * p.Q4 * Ca_n.^3; %[uM^-1]
            m_c         = (Ca_n ./ phi_N) .* dphi_N; %[-]
            CaM         = Ca_n ./ m_c; %[uM]
            tau_nk      = p.x_nk ^ 2 ./  (2 * p.D_cNO);
            p_NO_n      = p.NOswitch * (nNOS_act_n .* p.V_max_NO_n .* p.O2_n ./ (p.K_mO2_n + p.O2_n) * p.LArg_n / (p.K_mArg_n + p.LArg_n) ); %[uM/s] %*****
            c_NO_n      = p.k_O2_n * NO_n.^2 * p.O2_n; %[uM/s]
%             p_NO_n      = nNOS_act_n .* p.V_max_NO_n .* O2*1000 ./ (p.K_mO2_n + O2*1000) * p.LArg_n / (p.K_mArg_n + p.LArg_n); %[uM/s]
%             c_NO_n      = p.k_O2_n * NO_n.^2 .* O2*1000; %[uM/s]
            d_NO_n      = (NO_k - NO_n) ./ tau_nk; %[uM/s]
            
            %% Elshin neuron model
            %% Nernst potential for Na,K ions in soma and dendrite (Cl constant)
            E_Na_sa     = p.ph * log(Na_e ./ Na_sa);
            E_K_sa      = p.ph * log(K_e ./ K_sa);
            E_Na_d      = p.ph * log(Na_e ./ Na_d);
            E_K_d       = p.ph * log(K_e ./ K_d);
            
            %% Ion fluxes
            % Leak fluxes of Na,K,Cl in soma and dendrite using HH
            J_Naleak_sa = p.gNaleak_sa * (v_sa - E_Na_sa);
            J_Kleak_sa  = p.gKleak_sa * (v_sa - E_K_sa);
            J_Naleak_d  = p.gNaleak_d * (v_d - E_Na_d);
            J_Kleak_d   = p.gKleak_d * (v_d - E_K_d);

            % Na flux through NaP channel in soma using GHK
            m1alpha     = 1 ./ (6 * (1 + exp(-((0.143 * v_sa) + 5.67))));
            m1beta      = exp(-((0.143 * v_sa) + 5.67)) ./ (6 * (1 + exp(-((0.143 * v_sa) + 5.67))));
            h1alpha     = 5.12e-8 * exp(-((0.056 * v_sa) + 2.94));
            h1beta      = 1.6e-6 ./ (1 + exp(-(((0.2 * v_sa)) + 8)));
            J_NaP_sa    = (m1.^2 .* h1 .* p.gNaP_GHk * p.Farad .* v_sa .* (Na_sa - (exp(-v_sa / p.ph) .* Na_e))) ./ (p.ph * (1 - exp(-v_sa / p.ph)));

            % Na flux through NaT channel in soma using GHK
            m8alpha     = 0.32 * ((-v_sa - 51.9) ./ (exp(-(0.25 * v_sa + 12.975)) - 1));
            m8beta      = 0.28 * ((v_sa + 24.89) ./ (exp(0.2 * v_sa + 4.978) - 1));
            h6alpha     = 0.128 * exp(-(0.056 * v_sa + 2.94));
            h6beta      = 4 ./ (1 + exp(-(0.2 * v_sa + 6)));
            J_NaT_sa    = (m8.^3 .* h6 .* p.gNaT_GHk * p.Farad .* v_sa .* (Na_sa - (exp(-v_sa / p.ph) .* Na_e))) ./ (p.ph * (1 - exp(-v_sa / p.ph)));

            % K flux through KDR channel in soma using GHK
            m2alpha     = 0.016 * ((v_sa + 34.9) ./ (1 - exp(-((0.2 * v_sa) + 6.98))));
            m2beta      = 0.25 * exp(-((0.025 * v_sa) + 1.25));
            J_KDR_sa    =(m2.^2 .* p.gKDR_GHk * p.Farad .* v_sa .* (K_sa - (exp(-v_sa / p.ph) .* K_e))) ./ (p.ph * (1 - exp(-v_sa / p.ph)));

            % K flux through KA channel in soma using GHK
            m3alpha     = 0.02 * ((v_sa + 56.9) ./ (1 - exp(-((0.1 * v_sa) + 5.69))));
            m3beta      = 0.0175 * ((v_sa + 29.9) ./ (exp(((0.1 * v_sa) + 2.99)) - 1));
            h2alpha     = 0.016 * exp(-((0.056 * v_sa) + 4.61));
            h2beta      = 0.5 ./ (1 + exp(-((0.2 * v_sa) + 11.98)));
            J_KA_sa     = (m3.^2 .* h2 .* p.gKA_GHk * p.Farad .* v_sa .* (K_sa - (exp(-v_sa / p.ph) .* K_e))) ./ (p.ph * (1 - exp(-v_sa / p.ph)));

            % Na flux through NaP channel in dendrite using GHK
            m4alpha     = 1 ./ (6 * (1 + exp(-((0.143 * v_d) + 5.67))));
            m4beta      = exp(-((0.143 * v_d) + 5.67)) ./ (6 * (1 + exp(-((0.143 * v_d) + 5.67))));
            h3alpha     = 5.12e-8 * exp(-((0.056 * v_d) + 2.94));
            h3beta      = 1.6e-6 ./ (1 + exp(-(((0.2 * v_d)) + 8)));
            J_NaP_d     = (m4.^2 .* h3 .* p.gNaP_GHk * p.Farad .* v_d .* (Na_d - (exp(-v_d / p.ph) .* Na_e))) ./ (p.ph * (1 - exp(-v_d / p.ph)));

            % Na flux through NMDA channel in dendrite using GHK
            m5alpha     = 0.5 ./ (1 + exp((13.5 - K_e) / 1.42));
            m5beta      = 0.5 - m5alpha;
            h4alpha     = 1 ./ (2000 * (1 + exp((K_e - 6.75) / 0.71)));
            h4beta      = 1/2e3 - h4alpha;  %%% different from supplementary material!
            %h4beta      = 5e-5 - h4alpha;
            J_NMDA_d    = ( (m5 .* h4 .* p.gNMDA_GHk * p.Farad .* v_d .* (Na_d - (exp(-v_d / p.ph) .* Na_e))) ./ (p.ph * (1 - exp(-v_d / p.ph))) ) ./ (1 + 0.33 * p.Mg * exp(-(0.07 * v_d + 0.7)));

            % K flux through KDR channel in dendrite using GHK
            m6alpha     = 0.016 * ((v_d + 34.9) ./ (1 - exp(-((0.2 * v_d) + 6.98))));
            m6beta      = 0.25 * exp(-((0.025 * v_d) + 1.25));
            J_KDR_d     = (m6.^2 .* p.gKDR_GHk * p.Farad .* v_d .* (K_d - (exp(-v_d / p.ph) .* K_e))) ./ (p.ph * (1 - exp(-v_d / p.ph)));

            % K flux through KA channel in dendrite using GHK
            m7alpha     = 0.02 * ((v_d + 56.9) ./ (1 - exp(-((0.1 * v_d) + 5.69))));
            m7beta      = 0.0175 * ((v_d + 29.9) ./ (exp(((0.1 * v_d) + 2.99)) - 1));
            h5alpha     = 0.016 * exp(-((0.056 * v_d) + 4.61));
            h5beta      = 0.5 ./ (1 + exp(-((0.2 * v_d) + 11.98)));
            J_KA_d      = (m7.^2 .* h5 .* p.gKA_GHk * p.Farad .* v_d .* (K_d - (exp(-v_d / p.ph) .* K_e))) ./ (p.ph * (1 - exp(-v_d / p.ph)));
            
            % flux through the NaK-ATPase pump
            J_pump1_sa      = (1 + (p.K_init_e ./ K_e)).^(-2) .* (1 + (p.Na_init_sa ./ Na_sa)) .^ (-3);
            J_pump1init_sa  = (1 + (p.K_init_e / p.K_init_e)).^(-2) .* (1 + (p.Na_init_sa / p.Na_init_sa)).^(-3);
            J_pump1_d       = (1 + (p.K_init_e ./ K_e)).^(-2) .* (1 + (p.Na_init_d ./ Na_d)).^(-3);
            J_pump1init_d   = (1 + (p.K_init_e / p.K_init_e)).^(-2) .* (1 + (p.Na_init_d / p.Na_init_d)).^(-3);
            J_pump2         = 2 * (1 + p.O2_0 ./ (((1 - p.alph) * O2) + p.alph * p.O2_0)).^(-1);

            J_pump_sa   = p.Imax * J_pump1_sa .* J_pump2;
            J_pump_d    = p.Imax * J_pump1_d .* J_pump2;

            J_Napump_sa = 3 * J_pump_sa;
            J_Kpump_sa  = -2 * J_pump_sa;
            J_Napump_d  = 3 * J_pump_d;
            J_Kpump_d   = -2 * J_pump_d;  
            
            %% Total ion fluxes
            %Total ion fluxes in soma
            J_Na_tot_sa = J_NaP_sa + J_Naleak_sa + J_Napump_sa + p.NaTswitch *  J_NaT_sa;
            J_K_tot_sa  = J_KDR_sa + J_KA_sa + J_Kleak_sa + J_Kpump_sa;
            J_Cl_tot_sa = p.gClleak_sa * (v_sa- p.E_Cl_sa); % leak 
            
            % Total ion fluxes in dendrite
            J_Na_tot_d  = J_NaP_d + J_Naleak_d + J_Napump_d + J_NMDA_d;
            J_K_tot_d   = J_KDR_d + J_KA_d + J_Kleak_d + J_Kpump_d + J_NMDA_d;
            J_Cl_tot_d  = p.gClleak_d * (v_d - p.E_Cl_d); % leak
                        
            % Total ion fluxes in soma and dendrite
            J_tot_sa    = J_Na_tot_sa + J_K_tot_sa + J_Cl_tot_sa;
            J_tot_d     = J_Na_tot_d + J_K_tot_d + J_Cl_tot_d;
            
            % Fluxes for output
            J_K_NEtoECS_output = 1/(p.Farad * p.fe) * (((p.As * J_K_tot_sa) / p.Vs)  + ((p.Ad * J_K_tot_d) / p.Vd));
            
            J_K_NEtoSC_output =  1 / (p.Farad * p.fe) * ( ( (p.As * J_K_tot_sa) / p.Vs ) / p.VR_sn  +  ( (p.Ad * J_K_tot_d) / p.Vd ) / p.VR_dn ) * 1000;
            J_Na_NEtoSC_output = 1 / (p.Farad * p.fe) * ( ( (p.As * J_Na_tot_sa) / p.Vs ) / p.VR_sn  +  ( (p.Ad * J_Na_tot_d) / p.Vd ) / p.VR_dn ) * 1000;
            
            
            %% Tissue oxygen 
            P_02            = ( 2 * (1 + (p.O2_0 ./ ((1 - p.alph) * O2 + p.alph * p.O2_0))).^(-1) - (2 * (1 + (p.O2_0 / (p.alph * p.O2_0))).^(-1)) ) ./ ( (2 * (1 + (p.O2_0 / ((1 - p.alph) * p.O2_0 + p.alph * p.O2_0))).^(-1)) - (2 * (1 + (p.O2_0 / (p.alph * p.O2_0))).^(-1)) );  % P(02)
            CBF             = p.CBF_init * (R.^4 / p.LU_R_init^4);
            J_O2_vascular   = CBF .* ((p.O2_b - O2) ./ (p.O2_b - p.O2_0));
            J_O2_background = p.CBF_init * P_02 * (1 - p.gamm);
            J_O2_pump       = p.CBF_init * P_02 * p.gamm .* ((J_pump1_sa + J_pump1_d) ./ (J_pump1init_sa + J_pump1init_d));
            
            
        %% Conservation equations
            %% NO pathway
            du(idx.Ca_n, :)         = (I_Ca_tot / (2 * p.F * p.V_spine) - (p.k_ex * (Ca_n - p.Ca_rest))) / (1 + p.lambda_buf);  %[uM/s]
            du(idx.nNOS_act_n, :)   = p.V_maxNOS * CaM ./ (p.K_actNOS + CaM) - p.mu2_n * nNOS_act_n;                            %[uM/s]
            du(idx.NO_n, :)         = p_NO_n - c_NO_n + d_NO_n;                                                                 %[uM/s]
            
            %% Elshin neuron model
            
            % change in membrane potential
            du(idx.v_sa, :)     = 1/p.Cm * ( -J_tot_sa + 1 / (2 * p.Ra * p.dhod.^2) * (v_d - v_sa) + self.input_current(t) );
            du(idx.v_d, :)      = 1/p.Cm * (-J_tot_d + 1 / (2 * p.Ra * p.dhod.^2) * (v_sa - v_d));

            %change in concentration of Na,K,Cl in the soma
            du(idx.Na_sa, :)    = -p.As / (p.Farad * p.Vs) * J_Na_tot_sa + p.D_Na * (p.Vd + p.Vs) ./ (2 * p.dhod.^2 * p.Vs) * (Na_d - Na_sa);
            du(idx.K_sa, :)     = -p.As / (p.Farad * p.Vs) * J_K_tot_sa + p.D_K * (p.Vd + p.Vs) ./ (2 * p.dhod.^2 * p.Vs) * (K_d - K_sa);
            du(idx.Cl_sa, :)    = -p.As / (p.Farad * p.Vs) * J_Cl_tot_sa + p.D_Cl * (p.Vd + p.Vs) ./ (2 * p.dhod.^2 * p.Vs) * (Cl_d - Cl_sa);

            %change in concentration of Na,K,Cl in the dendrite
            du(idx.Na_d, :)     = -p.Ad / (p.Farad * p.Vd) * J_Na_tot_d + p.D_Na * (p.Vs + p.Vd) ./ (2 * p.dhod.^2 * p.Vd) * (Na_sa - Na_d);
            du(idx.K_d, :)      = -p.Ad / (p.Farad * p.Vd) * J_K_tot_d + p.D_K * (p.Vs + p.Vd) ./ (2 * p.dhod.^2 * p.Vd) * (K_sa - K_d);
            du(idx.Cl_d, :)     = -p.Ad / (p.Farad * p.Vd) * J_Cl_tot_d + p.D_Cl * (p.Vs + p.Vd) ./ (2 * p.dhod.^2 * p.Vd) * (Cl_sa - Cl_d);
            
            % change in potassium buffer for K+ in the extracellular space
%             J_buff              = p.Mu * K_e .* (p.B0 - Buff) ./ (1 + exp(-((K_e - 15) ./ 1.09))) - (p.Mu * Buff); % Elshin
            J_buff              = 8e-6 * K_e .* Buff .* exp((K_e - 5.5) ./ -1.09) - 8e-6 * (200 - Buff);   % Chang
            du(idx.Buff, :)     = J_buff;

            
            % change in concentration of Na,K,Cl in the extracellular space 
            du(idx.Na_e, :)     = 1/(p.Farad * p.fe) * (((p.As * J_Na_tot_sa) / p.Vs) + ((p.Ad * J_Na_tot_d) / p.Vd));
            du(idx.K_e, :)      = 1/(p.Farad * p.fe) * (((p.As * J_K_tot_sa) / p.Vs)  + ((p.Ad * J_K_tot_d) / p.Vd)) - p.buffSwitch*J_buff;
            du(idx.Cl_e, :)     = 1/(p.Farad * p.fe) * (((p.As * J_Cl_tot_sa) / p.Vs) + ((p.Ad * J_Cl_tot_d) / p.Vd));
            
            % change in tissue oxygen
            du(idx.O2, :)       = J_O2_vascular - J_O2_background - J_O2_pump;

            % Change in activation gating variables m
            du(idx.m1, :)       = 1000 * ((m1alpha .* (1 - m1)) - (m1beta .* m1));
            du(idx.m2, :)       = 1000 * ((m2alpha .* (1 - m2)) - (m2beta .* m2));
            du(idx.m3, :)       = 1000 * ((m3alpha .* (1 - m3)) - (m3beta .* m3));
            du(idx.m4, :)       = 1000 * ((m4alpha .* (1 - m4)) - (m4beta .* m4));
            du(idx.m5, :)       = 1000 * ((m5alpha .* (1 - m5)) - (m5beta .* m5));
            du(idx.m6, :)       = 1000 * ((m6alpha .* (1 - m6)) - (m6beta .* m6));
            du(idx.m7, :)       = 1000 * ((m7alpha .* (1 - m7)) - (m7beta .* m7));            
            du(idx.m8, :)       = 1000 * ((m8alpha .* (1 - m8)) - (m8beta .* m8));

            % Change in inactivation gating variables h
            du(idx.h1, :)       = 1000 * ((h1alpha .* (1 - h1)) - (h1beta .* h1));
            du(idx.h2, :)       = 1000 * ((h2alpha .* (1 - h2)) - (h2beta .* h2));
            du(idx.h3, :)       = 1000 * ((h3alpha .* (1 - h3)) - (h3beta .* h3));
            du(idx.h4, :)       = 1000 * ((h4alpha .* (1 - h4)) - (h4beta .* h4));
            du(idx.h5, :)       = 1000 * ((h5alpha .* (1 - h5)) - (h5beta .* h5));
            du(idx.h6, :)       = 1000 * ((h6alpha .* (1 - h6)) - (h6beta .* h6));

            
            du = bsxfun(@times, self.enabled, du);
            if nargout == 2
            	Uout = zeros(self.n_out, size(u, 2));
                Uout(self.idx_out.Glu, :) = self.input_Glu(t, v_sa);
                Uout(self.idx_out.current, :) = self.input_current(t);
            	Uout(self.idx_out.w_NR2A, :) = w_NR2A;
            	Uout(self.idx_out.w_NR2B, :) = w_NR2B;
            	Uout(self.idx_out.I_Ca, :) = I_Ca;
            	Uout(self.idx_out.I_Ca_tot, :) = I_Ca_tot;
            	Uout(self.idx_out.phi_N, :) = phi_N;
                Uout(self.idx_out.dphi_N, :) = dphi_N;
            	Uout(self.idx_out.m_c, :) = m_c;
            	Uout(self.idx_out.CaM, :) = CaM;
            	Uout(self.idx_out.tau_nk, :) = tau_nk;
            	Uout(self.idx_out.p_NO_n, :) = p_NO_n;
            	Uout(self.idx_out.c_NO_n, :) = c_NO_n;
            	Uout(self.idx_out.d_NO_n, :) = d_NO_n;
                Uout(self.idx_out.J_K_NEtoSC_output, :) = J_K_NEtoSC_output;
                Uout(self.idx_out.J_Na_NEtoSC_output, :) = J_Na_NEtoSC_output;
                Uout(self.idx_out.J_K_NEtoECS_output, :) = J_K_NEtoECS_output;
            	varargout = {Uout};
            end
        end
        
       function [J_K_NEtoSC, J_Na_NEtoSC, NO_n, v_sa] = shared(self, t, u, K_s, Na_s)    %shared variable
            p = self.params;
            idx = self.index;
            
            % Duplication of earlier equations - necessary but gross, fix
            % later7*10?5 S/cm2
            v_sa = u(idx.v_sa, :);      % membrane potential of soma/axon, mV
            v_d = u(idx.v_d, :);        % membrane potential of dendrite, mV
            K_sa = u(idx.K_sa, :);      % K+ concentration of soma/axon, mM
            Na_sa = u(idx.Na_sa, :);    % Na+ concentration of soma/axon, mM
            K_d = u(idx.K_d, :);        % K+ concentration of dendrite, mM
            Na_d = u(idx.Na_d, :);      % Na+ concentration of dendrite, mM
            K_e = u(idx.K_e, :);        % K+ concentration of ECS, mM
            Na_e = u(idx.Na_e, :);      % Na+ concentration of ECS, mM
            O2 = u(idx.O2, :);          % Tissue oxygen concentration, mM
            m1 = u(idx.m1, :);          % Activation gating variable, soma/axon NaP channel (Na+)
            m2 = u(idx.m2, :);          % Activation gating variable, soma/axon KDR channel (K+)
            m3 = u(idx.m3, :);          % Activation gating variable, soma/axon KA channel (K+)
            m4 = u(idx.m4, :);          % Activation gating variable, dendrite NaP channel (Na+)
            m5 = u(idx.m5, :);          % Activation gating variable, dendrite NMDA channel (Na+)
            m6 = u(idx.m6, :);          % Activation gating variable, dendrite KDR channel (K+)
            m7 = u(idx.m7, :);          % Activation gating variable, dendrite KA channel (K+)
            m8 = u(idx.m8, :);          % Activation gating variable, soma/axon NaT channel (Na+)
            h1 = u(idx.h1, :);          % Inactivation gating variable, soma/axon NaP channel (Na+)
            h2 = u(idx.h2, :);          % Inactivation gating variable, soma/axon KA channel (K+)
            h3 = u(idx.h3, :);          % Inactivation gating variable, dendrite NaP channel (Na+)
            h4 = u(idx.h4, :);          % Inactivation gating variable, dendrite NMDA channel (Na+)
            h5 = u(idx.h5, :);          % Inactivation gating variable, dendrite KA channel (K+)
            h6 = u(idx.h6, :);          % Inactivation gating variable, soma/axon NaT channel (Na+)
            Buff = u(idx.Buff, :);
            
            E_Na_sa     = p.ph * log(Na_e ./ Na_sa);
            E_K_sa      = p.ph * log(K_e ./ K_sa);
            E_Na_d      = p.ph * log(Na_e ./ Na_d);
            E_K_d       = p.ph * log(K_e ./ K_d);
            J_Naleak_sa = p.gNaleak_sa * (v_sa - E_Na_sa);
            J_Kleak_sa  = p.gKleak_sa * (v_sa - E_K_sa);
            J_Naleak_d  = p.gNaleak_d * (v_d - E_Na_d);
            J_Kleak_d   = p.gKleak_d * (v_d - E_K_d);
            J_NaP_sa    = (m1.^2 .* h1 .* p.gNaP_GHk * p.Farad .* v_sa .* (Na_sa - (exp(-v_sa / p.ph) .* Na_e))) ./ (p.ph * (1 - exp(-v_sa / p.ph)));
            J_NaT_sa    = (m8.^3 .* h6 .* p.gNaT_GHk * p.Farad .* v_sa .* (Na_sa - (exp(-v_sa / p.ph) .* Na_e))) ./ (p.ph * (1 - exp(-v_sa / p.ph)));
            J_KDR_sa    =(m2.^2 .* p.gKDR_GHk * p.Farad .* v_sa .* (K_sa - (exp(-v_sa / p.ph) .* K_e))) ./ (p.ph * (1 - exp(-v_sa / p.ph)));
            J_KA_sa     = (m3.^2 .* h2 .* p.gKA_GHk * p.Farad .* v_sa .* (K_sa - (exp(-v_sa / p.ph) .* K_e))) ./ (p.ph * (1 - exp(-v_sa / p.ph)));
            J_NaP_d     = (m4.^2 .* h3 .* p.gNaP_GHk * p.Farad .* v_d .* (Na_d - (exp(-v_d / p.ph) .* Na_e))) ./ (p.ph * (1 - exp(-v_d / p.ph)));
            J_NMDA_d    = ((m5 .* h4 .* p.gNMDA_GHk * p.Farad .* v_d .* (Na_d - (exp(-v_d / p.ph) .* Na_e))) ./ (p.ph * (1 - exp(-v_d / p.ph)))) ./ (1 + 0.33 * p.Mg * exp(-(0.07 * v_d + 0.7)));
            J_KDR_d     = (m6.^2 .* p.gKDR_GHk * p.Farad .* v_d .* (K_d - (exp(-v_d / p.ph) .* K_e))) ./ (p.ph * (1 - exp(-v_d / p.ph)));
            J_KA_d      = (m7.^2 .* h5 .* p.gKA_GHk * p.Farad .* v_d .* (K_d - (exp(-v_d / p.ph) .* K_e))) ./ (p.ph * (1 - exp(-v_d / p.ph)));
            J_pump1_sa  = (1 + (p.K_init_e ./ K_e)).^(-2) .* (1 + (p.Na_init_sa ./ Na_sa)) .^ (-3);
            J_pump1_d   = (1 + (p.K_init_e ./ K_e)).^(-2) .* (1 + (p.Na_init_d ./ Na_d)).^(-3);
            J_pump2     = 2 * (1 + p.O2_0 ./ (((1 - p.alph) * O2) + p.alph * p.O2_0)).^(-1);
            J_pump_sa   = p.Imax * J_pump1_sa .* J_pump2;
            J_pump_d    = p.Imax * J_pump1_d .* J_pump2;
            J_Napump_sa = 3 * J_pump_sa;
            J_Kpump_sa  = -2 * J_pump_sa;
            J_Napump_d  = 3 * J_pump_d;
            J_Kpump_d   = -2 * J_pump_d;             
            J_Na_tot_d  = J_NaP_d + J_Naleak_d + J_Napump_d + J_NMDA_d;
            J_K_tot_d   = J_KDR_d + J_KA_d + J_Kleak_d + J_Kpump_d + J_NMDA_d;
            J_Na_tot_sa = J_NaP_sa + J_Naleak_sa + J_Napump_sa + p.NaTswitch * J_NaT_sa;
            J_K_tot_sa  = J_KDR_sa + J_KA_sa + J_Kleak_sa + J_Kpump_sa;
            
%             J_buff              = p.Mu * K_e .* (p.B0 - Buff) ./ (1 + exp(-((K_e - 15) ./ 1.09))) - (p.Mu * Buff); % Elshin
            J_buff              = 8e-6 * K_e .* Buff .* exp((K_e - 5.5) ./ -1.09) - 8e-6 * (200 - Buff);   % Chang
            
            %Fluxes from neuron to SC
%             J_K_NEtoSC =  1 / (p.Farad * p.fe) * ( ( (p.As * J_K_tot_sa) / p.Vs ) / p.VR_sn  +  ( (p.Ad * J_K_tot_d) / p.Vd ) / p.VR_dn ) * 1000;
%             J_Na_NEtoSC = 1 / (p.Farad * p.fe) * ( ( (p.As * J_Na_tot_sa) / p.Vs ) / p.VR_sn  +  ( (p.Ad * J_Na_tot_d) / p.Vd ) / p.VR_dn ) * 1000;
              J_K_NEtoSC =  1 / (p.Farad * p.fe) * ( ( (p.As* J_K_tot_sa) / p.Vs )   +  ( (p.Ad * J_K_tot_d) / p.Vd  )  ) * 1000 - p.buffSwitch*J_buff*1000;
              J_Na_NEtoSC = 1 / (p.Farad * p.fe) * ( ( (p.As* J_Na_tot_sa) / p.Vs  )  +  ( (p.Ad * J_Na_tot_d) / p.Vd  ) ) * 1000;            

            NO_n = u(idx.NO_n, :);
            
       end
        % This is basically a scaled version of input_rho(t) in Astrocyte
        function Glu = input_Glu(self, t, v_sa) 
            p = self.params;
%             Glu = p.GluSwitch * (p.Glu_max - p.Glu_min) * ( ...
%             0.5 * tanh((t - p.t_0_Glu) / p.theta_L_Glu) - ...
%             0.5 * tanh((t - p.t_2_Glu) / p.theta_R_Glu)) + p.Glu_min;
            vsaMin = -70; vsaMax = -20;
            if (v_sa > vsaMax)
                Glu = p.GluSwitch * 1846;
            elseif (v_sa < vsaMin)
                Glu = p.GluSwitch * 0;
            else
                Glu = p.GluSwitch * ( 1846/(vsaMax - vsaMin) ) * (v_sa - vsaMin);
            end

        end
        
        % Current input to neuron
        function current = input_current(self, t) 
            p = self.params;
            t_0 = p.t_0;
            t_2 = p.t_2;        
            
            if p.currentType == 1
                current = 0.03*gaussmf(t, [6, t_0+15]);
            elseif p.currentType == 2
                current = 0.03 * (0.5 * tanh((t - t_0)/5) - ...
                0.5 * tanh((t - t_2)/5));
            end
        end        
        
        function names = varnames(self)
            names = [fieldnames(self.index); fieldnames(self.idx_out)];
        end
    end
end    
        
function idx = indices()    %for state variables
    idx.Ca_n = 1;                  
    idx.nNOS_act_n = 2;
    idx.NO_n = 3;
    idx.v_sa = 4;
    idx.v_d = 5;
    idx.K_sa = 6;
    idx.Na_sa = 7;
    idx.Cl_sa = 8;
    idx.K_d = 9;
    idx.Na_d = 10;
    idx.Cl_d = 11;
    idx.K_e = 12;
    idx.Na_e = 13;
    idx.Cl_e = 14;
    idx.Buff = 15;
    idx.O2 = 16;
    idx.m1 = 17;
    idx.m2 = 18;
    idx.m3 = 19;
    idx.m4 = 20;
    idx.m5 = 21;
    idx.m6= 22;
    idx.m7 = 23;
    idx.m8 = 24;
    idx.h1 = 25;
    idx.h2 = 26;
    idx.h3 = 27;
    idx.h4 = 28;
    idx.h5 = 29;
    idx.h6 = 30;
            
end        

function [idx, n] = output_indices()    %for variables in nargout loop
    idx.ft = 1;
    idx.J_Na_n = 2;
    idx.Glu = 3; 
    idx.w_NR2A = 4; 
    idx.w_NR2B = 5; 
    idx.I_Ca = 6; 
    idx.I_Ca_tot = 7; 
    idx.phi_N = 8; 
    idx.dphi_N = 9; 
    idx.m_c = 10; 
    idx.CaM = 11; 
    idx.tau_nk = 12; 
    idx.p_NO_n = 13;    
    idx.c_NO_n = 14; 
    idx.d_NO_n = 15;    
    idx.current = 16;
    idx.J_K_NEtoSC_output = 18;
    idx.J_Na_NEtoSC_output = 19;
    idx.J_K_NEtoECS_output = 20;
    n = numel(fieldnames(idx));
end
        
function params = parse_inputs(varargin)
    parser = inputParser();
    
    parser.addParameter('GluSwitch', 1); 
    parser.addParameter('NaTswitch', 1); 
    parser.addParameter('buffSwitch', 1); 
    parser.addParameter('NOswitch', 1); 
    
    parser.addParameter('currentType', 1); 
    
    % global constants
    parser.addParameter('F', 9.65e4); %C mol^-1; Faraday's constant
    parser.addParameter('R_gas', 8.315); %J mol^-1 K^-1; Gas constant
    parser.addParameter('T', 310); % K; Temperature
    
    % input 
    parser.addParameter('startpulse', 200);    
    parser.addParameter('lengthpulse', 200);
    parser.addParameter('lengtht1', 10);
    parser.addParameter('F_input', 2.5);        %s  
    parser.addParameter('alpha', 2);
    parser.addParameter('beta', 5);
    parser.addParameter('delta_t', 10);         %s
    parser.addParameter('k_C', 7.35e-5);        %uM m s^-1
    parser.addParameter('Glu_max', 1846);       % microM (one vesicle, Santucci2008)
    parser.addParameter('Glu_min', 0);          % microM
    parser.addParameter('theta_L_Glu', 1);      % slope of Glu input 
    parser.addParameter('theta_R_Glu', 1);      % slope of Glu input 

    % NO pathway
    parser.addParameter('K_mA', 650);           % [uM] - fit to Santucci2008
    parser.addParameter('K_mB', 2800);          % [uM] - fit to Santucci2008

    parser.addParameter('v_n', -0.04);          % [V] ; the neuronal membrane potential , assumed to be approx constant in this model ********
    parser.addParameter('G_M', 46000);          % [fS]! was 46 pS! ; the conductance of the NMDA channel to Ca2+ compaired  
    parser.addParameter('P_Ca_P_M', 3.6);       % [-] ; the relative conductance of the NMDA channel to Ca2+ compared to monovalent ions
    parser.addParameter('Ca_ex', 2e3);          % [microM] ; the external calcium concentration (in Comerford+David2008: 1.5 mM!)
    parser.addParameter('M', 1.3e5);            % [microM] ; the concentration of monovalent ions in the neuron

    parser.addParameter('n_NR2A', 0.63);        % [-] ; average number of NR2A NMDA receptors per synapse (Santucci2008)
    parser.addParameter('n_NR2B', 11);          % [-] ; average number of NR2B NMDA receptors per synapse (Santucci2008)
    parser.addParameter('V_max_NO_n', 4.22);    % [s^-1] ; maximum catalytic rate of NO production (Chen2006) - obtained from fig 6 & equ 17 & 18
    parser.addParameter('O2_n', 200);           % [uM] ; tissue O2 concentration in the neuron (M.E.) *************
    parser.addParameter('K_mO2_n', 243);        % [uM] ; Chen2006
    parser.addParameter('LArg_n', 100);         % [uM] ; 
    parser.addParameter('K_mArg_n', 1.5);       % [uM] ; 
    
    parser.addParameter('k_O2_n', 9.6e-6);      % [uM^-2 s^-1] ; % (Kavdia2002)
    parser.addParameter('x_nk', 25);            % [um] ;  (M.E.)
    parser.addParameter('D_cNO', 3300);         % [um^2 s^-1] ; Diffusion coefficient NO (Malinski1993)
    
    parser.addParameter('V_spine', 8e-8);       % [nL] ; volume of the neuronal dendritic spine Santucci2008
    parser.addParameter('k_ex', 1600);          % [s^-1] ; decay rate constant of internal calcium concentration Santucci2008
    parser.addParameter('Ca_rest', 0.1);        % [uM] ; resting calcium concentration (in Comerford+David2008: 2.830 mM; in Santucci2008P: 0.1 \muM)
    parser.addParameter('lambda_buf', 20);      % [-] ; buffer capacity Santucci2008

    parser.addParameter('V_maxNOS', 25e-3);     % [] ; M.E.
    parser.addParameter('K_actNOS', 9.27e-2);   % [uM] ; 
    parser.addParameter('mu2_n', 0.0167);       % [s^-1] ; rate constant at which the nNOS is deactivated Comerford2008
    
    parser.addParameter('Q1', 1.9e5);           % [uM^-1]
    parser.addParameter('Q2', 2.1e5);           % [uM^-1]
    parser.addParameter('Q3', 0.4e5);           % [uM^-1]
    parser.addParameter('Q4', 0.26e5);          % [uM^-1]

    % Elshin neuron model
    parser.addParameter('E_Cl_sa', -70);        % Nernst potential for Cl- in soma      [mV]
    parser.addParameter('E_Cl_d', -70);         % Nernst potential for Cl- in dendrite   [mV]
    parser.addParameter('Ra', 1.83e5);          % Input resistance of dendritic tree [ohms]
    parser.addParameter('dhod', 4.5e-2);        % Half length of dendrite [cm]
    parser.addParameter('As', 1.586e-5);        % Surface area of soma [cm2]
    parser.addParameter('Ad', 2.6732e-4);       % Surface area of dendrite [cm2]
    parser.addParameter('Vs', 2.16e-9);         % Volume of soma [cm3]
    parser.addParameter('Vd', 5.614e-9);        % Volume of dendrite [cm3]
    parser.addParameter('VR_sn', 0.2778);       % Percentage volume of soma in neuron Vs / (Vs + Vd)
    parser.addParameter('VR_dn', 0.7222);       % Percentage volume of dendrite in neuron Vd / (Vs + Vd)    
    parser.addParameter('fe', 0.15);            % ECS to neuron ratio
    parser.addParameter('Cm', 7.5e-7);          % Membrance capacitance [s/ohms cm2]    *******
%     parser.addParameter('Cm', 0.75);          % Membrance capacitance [uF / cm2]    *******
    parser.addParameter('Farad', 96.485);           % Faraday constant [C / mmol]
    parser.addParameter('ph', 26.6995);         % RT/F
    parser.addParameter('Mu', 8e-4);            % [m/s]
    parser.addParameter('B0', 500);             % Effective total buffer concentration [mM]
    
    parser.addParameter('gNaP_GHk', 2e-6);
    parser.addParameter('gKDR_GHk', 10e-5);
    parser.addParameter('gKA_GHk', 1e-5);
    parser.addParameter('gNMDA_GHk', 1e-5);   
    parser.addParameter('gNaT_GHk', 10e-5); 
    
%     parser.addParameter('gNaleak_sa', 4.0678e-5);     % Elshin's new values ******
%     parser.addParameter('gKleak_sa', 1.4148e-4);
%     parser.addParameter('gClleak_sa', 10*4.0678e-5);
%     parser.addParameter('gNaleak_d', 4.135e-5);   
%     parser.addParameter('gKleak_d', 1.4148e-4); 
%     parser.addParameter('gClleak_d', 10*4.135e-5);  
%     parser.addParameter('Imax', 0.013*4);             
    parser.addParameter('Imax', 1.48e-3);               % Chang values ******
    parser.addParameter('gNaleak_sa', 9.5999e-6 );
    parser.addParameter('gKleak_sa', 3.4564e-5);
    parser.addParameter('gClleak_sa', 10*9.5999e-6);
    parser.addParameter('gNaleak_d', 1.0187e-5 );   
    parser.addParameter('gKleak_d', 3.4564e-5); 
    parser.addParameter('gClleak_d', 10*1.0187e-5);   
    
    parser.addParameter('O2_0', 2e-2);          % [mM]
    parser.addParameter('alph',  0.05);         % percentage of ATP production independent of O2
    parser.addParameter('D_Na', 1.33e-5);       % [cm2 / s]
    parser.addParameter('D_K', 1.96e-5);        % [cm2 / s]
    parser.addParameter('D_Cl', 2.03e-5);       % [cm2 / s]
    
    parser.addParameter('K_init_e', 2.9); 
    parser.addParameter('Na_init_sa', 10); 
    parser.addParameter('Na_init_d', 10); 
    
    parser.addParameter('LU_R_init', 1.9341e-5); 
    parser.addParameter('CBF_init', 3.2e-2); 
    parser.addParameter('O2_b', 4e-2);   
    parser.addParameter('gamm', 0.10);  
    parser.addParameter('Mg', 1.2);  
    
    parser.parse(varargin{:})
    params = parser.Results;
    params.t_0 = params.startpulse;
    params.t_1 = params.t_0 + params.lengtht1;
    params.t_2 = params.t_0 + params.lengthpulse;
    params.t_3 = params.t_1 + params.lengthpulse;
    params.gab = factorial(params.alpha + params.beta - 1);
    params.ga = factorial(params.alpha - 1);
    params.gb = factorial(params.beta - 1);
    params.t_0_Glu = params.t_0;
    params.t_2_Glu = params.t_2;
    params.startpulse_Glu = params.startpulse;
    params.lengthpulse_Glu = params.lengthpulse;
end

function u0 = initial_conditions(idx)
    u0 = zeros(length(fieldnames(idx)), 1);
    
    u0(idx.Ca_n) = 0.0001;
    u0(idx.nNOS_act_n) = 0.3;
    u0(idx.NO_n) = 0.1;
    u0(idx.v_sa) = -70;
    u0(idx.v_d) = -70;
    u0(idx.K_sa) = 133.5;
    u0(idx.Na_sa) = 9.9854;
    u0(idx.Cl_sa) = 10.464;
    u0(idx.K_d) = 133.5;
    u0(idx.Na_d) = 9.9853;
    u0(idx.Cl_d) = 10.464;
    u0(idx.K_e) = 3.5006;
    u0(idx.Na_e) = 139.76;
    u0(idx.Cl_e) = 144.1;
    u0(idx.Buff) = 500;
    u0(idx.O2) = 0.022715;
    u0(idx.m1) = 0.012869;
    u0(idx.m2) = 0.0012175;
    u0(idx.m3) = 0.1193;
    u0(idx.m4) = 0.012869;
    u0(idx.m5) = 0.00087377;
    u0(idx.m6) = 0.0012175;
    u0(idx.m7) = 0.1193;
    u0(idx.m8) = 0.005;
    u0(idx.h1) = 0.9718;
    u0(idx.h2) = 0.12053;
    u0(idx.h3) = 0.9718;
    u0(idx.h4) = 0.99005;
    u0(idx.h5) = 0.12053;
    u0(idx.h6) = 0.9961;

end
