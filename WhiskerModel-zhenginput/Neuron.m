classdef Neuron < handle
    properties
        input_data
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
        function [du, varargout] = rhs(self, t, u, NO_k, R)
            t = t(:).';
            p = self.params;
            idx = self.index;
            
    %% State variables
            E_t =   u(idx.E_t, :);              % [-]                       Proportion of excitatory cells firing
            I_t =   u(idx.I_t, :);              % [-]                       Proportion of inhibitory cells firing
            K_e =   u(idx.K_e, :);              % [mM]                      K+ concentration of ECS
            Na_sa = u(idx.Na_sa, :);            % [mM]                      Na+ concentration of soma/axon
            Na_d =  u(idx.Na_d, :);             % [mM]                      Na+ concentration of dendrite 
            O2 =    u(idx.O2, :);               % [mM]                      Tissue oxygen concentration
            CBV =   u(idx.CBV, :);              % [-]                       BOLD CBV
            HbR =   u(idx.HbR, :);              % [-]                       BOLD deoxyhemoglobin
            Ca_n =  u(idx.Ca_n, :);             % [uM]                      Calcium concentration     
            nNOS_act_n = u(idx.nNOS_act_n, :);  % [uM]                      NO active concentration
            NO_n =  u(idx.NO_n, :);             % [uM]                      NO concentration
            
    %% Algebraic Expressions
            % E and I cells
            arg_e = p.c1 .* E_t - p.c2 .* I_t + p.w_te .* self.T(t);
                                                % [-]                       Excitatory cells response function input
            arg_i = p.c3 .* E_t - p.c4 .* I_t + p.w_ti .* self.T(t); 
                                                % [-]                       Inhibitory cells response function input
            
            % Oxygen 
            % O2switch=0 ATP is plentiful, O2switch=1 ATP is limited (oxygen-limited regime)
            O2_p            = p.O2_0 * (1 - p.O2switch) + O2 * p.O2switch;
                                                % [mM]                      Oxygen concentration dependent on ATP
            J_pump2         = 2 * (1 + p.O2_0 ./ (((1 - p.alpha_O2) * O2_p)...
                 + p.alpha_O2 * p.O2_0)).^(-1); % [mM s^-1]
            
            J_pump2_0       = 0.0952;           % [mM s^-1]                      
            J_pump2_O2_0    = 1;                % [mM s^-1]                         
            P_02            = (J_pump2 - J_pump2_0 ) ./ ( J_pump2_O2_0 -...
                J_pump2_0);                     % [-]                       Normalised pump rate of oxygen            
            
            CBF             = p.CBF_init * (R.^4 / p.R_init^4);
                                                % [-]                       Normalised cerebral bloof flow
            J_O2_vascular   = CBF .* ((p.O2_b - O2) ./ (p.O2_b - p.O2_0));
                                                % [mM s^-1]                 Vascular supply of oxygen
            J_O2_background = p.CBF_init * P_02 * (1 - p.gamma_O2);
                                                % [mM s^-1]                 Background oxygen consumption
                                                
            J_pump1_sa      = (1 + (p.K_init_e ./ K_e)).^(-2) .* (1 + ...
                (p.Na_init_sa ./ Na_sa)) .^ (-3); % [mM s^-1]               Oxygen consumption in soma/axon pump
            J_pump1init_sa  = 0.0312;           % [mM s^-1]                 Initial consumption in soma/axon pump    
            J_pump1_d       = (1 + (p.K_init_e ./ K_e)).^(-2) .* (1 + ...
                (p.Na_init_d ./ Na_d)).^(-3);   % [mM s^-1]                 Oxygen consumption in dendrite pump
            J_pump1init_d   = 0.0312;           % [mM s^-1]                 Initial consumption in dendrite pump            
            J_O2_pump       = p.CBF_init * P_02 * p.gamma_O2 .* ((J_pump1_sa ...
                + J_pump1_d) ./ (J_pump1init_sa + J_pump1init_d));
                                                % [mM s^-1]                 Total oxygen consumption of the NaK pumps
 
                                                
            % BOLD
            f_out           = CBV.^(1/p.d) + p.tau_TAT * (1/(p.tau_MTT + ...
                p.tau_TAT) .* ( CBF/p.CBF_init  - CBV.^(1/p.d) ));
                                                % [-]                       Buxton balloon model outflow
            CMRO2           = J_O2_background + J_O2_pump;
                                                % [mM s^-1]                 Metabolic rate of oxygen in the neurons
            CMRO2_init      = p.CBF_init * P_02;% [mM s^-1]                 Initial metabolic rate of oxygen
            
            
            % NO pathway
            % Glutamate input: vesicle released when the extracellular K+ is over 5.5 mM (Ke_switch)
            Glu = self.shared(t, u);            % [mM]

            w_NR2A = Glu ./ (p.K_mA + Glu);     % [-] 
            w_NR2B = Glu ./ (p.K_mB + Glu);     % [-]
          
            I_Ca = (-4 * p.v_n * p.G_M * p.P_Ca_P_M * (p.Ca_ex / p.M)) / ...
                (1 + exp(-0.08 * (p.v_n + 20))) * (exp(2 * p.v_n / p.ph)) / (1 - exp(2 * p.v_n / p.ph)); 
                                                % [fA]
            I_Ca_tot = I_Ca .* (p.n_NR2A * w_NR2A + p.n_NR2B * w_NR2B); 
                                                % [fA]
            
            CaM = Ca_n / p.m_c;                 % [uM]
            tau_nk = p.x_nk ^ 2 ./  (2 * p.D_cNO);
            
            p_NO_n = p.NOswitch * ( nNOS_act_n * p.V_max_NO_n * p.O2_n / ...
                (p.K_mO2_n + p.O2_n) * p.LArg_n / (p.K_mArg_n + p.LArg_n) );
            c_NO_n = p.k_O2_n * NO_n.^2 * p.O2_n;                           
            d_NO_n = (NO_k - NO_n) ./ tau_nk;                              
            
            % NEtoSC; the neuronal activation in the Dormanns simulation
            dKedt = -p.beta_K_e * (K_e - p.K_eBase) + ...
                (p.alpha_K_e .* p.beta_K_e) .* ((E_t) ./ p.E_relative);
            J_K_NEtoSC = p.k_syn * dKedt;                                 
            J_Na_NEtoSC = -p.k_syn * dKedt;
            
      %% Derivative equations                                                        
            
            %% Wilson & Cowan model
            % Excitatory and inhibitory cells (Input of this model)
            du(idx.E_t, :) = 1 ./ p.tau_e .* (-E_t + self.s_e(arg_e));
                                                % [ms^-1]                   Change in excitatory cells firing
            du(idx.I_t, :) = 1 ./ p.tau_i .* (-I_t + self.s_i(arg_i));
                                                % [ms^-1]                   Change in inhibitory cells firing
    
            % K_e, Na_sa and Na_d
            du(idx.K_e, :) = - p.beta_K_e * (K_e - p.K_eBase) + ...         
                p.alpha_K_e * p.beta_K_e .* ((E_t) ./ p.E_relative);
                                                % [mM ms^-1]                Change in extracellular potassium concentration
            du(idx.Na_sa, :) = - p.beta_Na_sa * (Na_sa - p.Na_saBase) + ...
                p.alpha_Na_sa * p.beta_Na_sa .* ((E_t) ./ p.E_relative);
                                                % [mM ms^-1]                Change in soma/axon sodium concentration
            du(idx.Na_d, :) = - p.beta_Na_d * (Na_d - p.Na_dBase) + ...
                p.alpha_Na_d * p.beta_Na_d .* ((E_t) ./ p.E_relative);
                                                % [mM ms^-1]                Change in dendrite sodium concentration
    
            % O2
            du(idx.O2, :) = J_O2_vascular - J_O2_background - J_O2_pump;
                                                % [uM ms^-1]                Change in oxygen concentration
            
            % BOLD
            du(idx.CBV, :) = 1/(p.tau_MTT + p.tau_TAT) .* ( CBF/p.CBF_init  - CBV.^(1/p.d) );
                                                % [ms^-1]                   Change in normalized cerebral blood volume
            du(idx.HbR, :) = 1/p.tau_MTT * ( CMRO2./CMRO2_init - HbR./CBV .* f_out );
                                                % [ms^-1]                   Change in normalized deoxyhemoglobin
    
            % NO
            du(idx.Ca_n, :) = (I_Ca_tot / (2 * p.Farad * p.V_spine) ...
                - (p.k_ex * (Ca_n - p.Ca_rest))) / (1 + p.lambda_buf);
                                                % [uM ms^-1]                Change in calcium concentration
            du(idx.nNOS_act_n, :) = p.V_maxNOS * CaM ./ (p.K_actNOS + CaM)...
                - p.mu2_n * nNOS_act_n;         % [uM ms^-1]                Change in active NO concentration
            du(idx.NO_n, :) = p_NO_n - c_NO_n + d_NO_n;
                                                % [uM ms^-1]                Change in NO concentration
            
            
            du = bsxfun(@times, self.enabled, du);
            if nargout == 2
            	Uout = zeros(self.n_out, size(u, 2));
                Uout(self.idx_out.Glu, :) = Glu;
                Uout(self.idx_out.dKedt, :) = dKedt;
                Uout(self.idx_out.J_K_NEtoSC, :) = J_K_NEtoSC;
                Uout(self.idx_out.J_Na_NEtoSC, :) = J_Na_NEtoSC;
                Uout(self.idx_out.CBF, :) = CBF;
                Uout(self.idx_out.CMRO2, :) = CMRO2;               
            	varargout = {Uout};
            end
        end

        function [Glu, J_K_NEtoSC, J_Na_NEtoSC, NO_n, O2] = shared(self, ~, u)
            p = self.params;
            idx = self.index;
            E_t = u(idx.E_t, :);
            I_t = u(idx.I_t, :);
            K_e = u(idx.K_e, :);
            O2 = u(idx.O2, :);
            NO_n = u(idx.NO_n, :);
            
            %% Glutamate function dependent on K_e (ECS K+)
            Glu = p.GluSwitch * 0.5 * p.Glu_max * ( 1 + tanh( (K_e - p.Ke_switch)...
                / p.Glu_slope) );               % [mM]                      Glutamate concentration            
            
            dKedt = -p.beta_K_e * (K_e - p.K_eBase) + ...
                (p.alpha_K_e .* p.beta_K_e) .* ((E_t) ./ p.E_relative);    % K_e derivative as neuron activation in Ostby
            J_K_NEtoSC = p.k_syn * dKedt;                                   % Input current: K+ flux into SC
            J_Na_NEtoSC = -p.k_syn * dKedt;                                 % Input current: Na+ flux out of SC
        end       
        
        function sigmoid_e = s_e(self,x)                                    % Response function for excitatory neurons
            p = self.params;
            epsilon_e = 5.12;
            sigmoid_e = epsilon_e./(1+exp(-(x-p.theta_e)./p.a_e));
        end
        function sigmoid_i = s_i(self,x)                                    % Response function for inhibitory neurons
            p = self.params;
            epsilon_i = 11.61;
            sigmoid_i = epsilon_i./(1+exp(-(x-p.theta_i)./p.a_i));
        end
        function input_e = T(self,t)                                        % External thalamic input given to the neurons
            p = self.params;
            if p.Whiskerpad == 0
                pulse_width_P = p.startpulse + p.lengthpulse;
                input_e = 1.0 .* (heaviside(t - p.startpulse) - heaviside(t - pulse_width_P));
            elseif p.Whiskerpad == 1 || p.Whiskerpad == 2
                neural_data = self.input_data;
                index_t = round(t + 1);
                input_e = neural_data(index_t);
            end
        end
        function names = varnames(self)
            names = [fieldnames(self.index); fieldnames(self.idx_out)];
        end
    end
end    
        
function idx = indices()
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
end        

function [idx, n] = output_indices()
    idx.Glu = 1;
    idx.dKedt = 2;
    idx.J_K_NEtoSC = 3;
    idx.J_Na_NEtoSC = 4;
    idx.CBF = 5;
    idx.CMRO2 = 6;
    n = numel(fieldnames(idx));     
end
      
function params = parse_inputs(varargin)
    parser = inputParser();
    
    % On and off switches for Glu, NO and O2
    parser.addParameter('GluSwitch', 1);  
    parser.addParameter('NOswitch', 1); 
    parser.addParameter('O2switch', 1);
    
    % Wilson and Cowan model parameters
    % K_e, Na_sa, Na_d parameters
    parser.addParameter('K_eBase', 3.567);      % [mM]
    parser.addParameter('Na_saBase', 9.237);    % [mM]
    parser.addParameter('Na_dBase', 9.318);     % [mM]
    
    parser.addParameter('alpha_K_e', 3.5);      % [-]                
    parser.addParameter('beta_K_e', 10e-3);     % [ms^-1]
    parser.addParameter('alpha_Na_sa', 4.23);   % [-]
    parser.addParameter('beta_Na_sa', 0.65e-3); % [ms^-1]
    parser.addParameter('alpha_Na_d', -2.12);   % [-]                           
    parser.addParameter('beta_Na_d', 0.7e-3);   % [ms^-1]
    
    % E and I parameters
    parser.addParameter('E_relative', 4.8923); % [-]                       Maximum value for E after external input
    parser.addParameter('c1', 42);              % [-]
    parser.addParameter('c2', 24.6);            % [-]
    parser.addParameter('c3', 42);              % [-]
    parser.addParameter('c4', 18);              % [-]
    parser.addParameter('w_te', 53.43);         % [-]
    parser.addParameter('w_ti', 68.4);          % [-]
    parser.addParameter('tau_e', 5);            % [ms]
    parser.addParameter('tau_i', 15);           % [ms]
    
    parser.addParameter('a_e', 4.16);           % [-]
    parser.addParameter('theta_e', 15);         % [-]
    parser.addParameter('a_i', 3.94);           % [-]
    parser.addParameter('theta_i', 15);         % [-]
    
    parser.addParameter('Pinput', 1.0);         % [-]
    parser.addParameter('Qinput', 0.0);         % [-]
    
    % Glutamate parameters
    parser.addParameter('Glu_max', 1846);       % [uM] (one vesicle, Santucci2008)
    parser.addParameter('Glu_slope', 0.1);      % [mM] Slope of sigmoidal (M.E.)
    parser.addParameter('Ke_switch', 5.5);      % [mM] Threshold past which 
    % glutamate vesicle is released (M.E.)   
    
    % Elshin model constants
    parser.addParameter('k_syn', 11.5);         % [-] 
    % Coupling scaling factor, values can range between 1.06 to 14.95 
    % according to estimation from experimental data of Ventura and Harris, 
    % Maximum dilation at 9.5
       
    % BOLD constants
    parser.addParameter('tau_MTT', 3e3);        % [ms]
    parser.addParameter('tau_TAT', 20e3);       % [ms]
    parser.addParameter('d', 0.4);              % [-]
    parser.addParameter('a_1', 3.4);            % [-]   
    parser.addParameter('a_2', 1);              % [-]
    parser.addParameter('V_0', 0.03);           % [-]
    parser.addParameter('E_0', 0.4);            % [-]
        
    % Elshin neuron model
    parser.addParameter('Farad', 96.485);       % [C / mmol]
    parser.addParameter('ph', 26.6995);         % [RT/F]

    % Oxygen model

    parser.addParameter('O2_0', 0.02);          % [mM]
    parser.addParameter('alpha_O2',  0.05);     % [-] Percentage of ATP 
    
    % production independent of O2
    parser.addParameter('K_init_e', 2.9);       % [mM]
    parser.addParameter('Na_init_sa', 10);      % [mM]
    parser.addParameter('Na_init_d', 10);       % [mM]
    
    parser.addParameter('R_init', 20);          % [um]
    parser.addParameter('CBF_init', 3.2e-2);    % [-]
    parser.addParameter('O2_b', 4e-2);          % [mM]
    parser.addParameter('gamma_O2', 0.10);      % [-]

    % NO pathway 
    parser.addParameter('m_c', 4);              % [-] 
    % Number of Ca2+ bound per calmodulin (approximated as parameter, 
    % originally an algebraic variable that changed from 3.999 to 4)
    parser.addParameter('K_mA', 650);           % [uM] - fit to Santucci2008
    parser.addParameter('K_mB', 2800);          % [uM] - fit to Santucci2008
    parser.addParameter('v_n', -40);            % [mV] ; the neuronal membrane 
    % potential , assumed to be approx constant in this model
    parser.addParameter('G_M', 46e9);           % [pS] = [kg^-1*m^-2*ms^3*A^2] --> converted to ms
    % channel to Ca2+ compaired  
    parser.addParameter('P_Ca_P_M', 3.6);       % [-] ; the relative conductance 
    % of the NMDA channel to Ca2+ compared to monovalent ions
    parser.addParameter('Ca_ex', 2e3);          % [uM] ; the external calcium 
    % concentration (in Comerford+David2008: 1.5 mM!)
    parser.addParameter('M', 1.3e5);            % [uM] ; the concentration of 
    % monovalent ions in the neuron
    parser.addParameter('n_NR2A', 0.63);        % [-] ; average number of 
    % NR2A NMDA receptors per synapse (Santucci2008)
    parser.addParameter('n_NR2B', 11);          % [-] ; average number of 
    % NR2B NMDA receptors per synapse (Santucci2008)
    parser.addParameter('V_max_NO_n', 4.22e-3); % [ms^-1] ; maximum catalytic
    % rate of NO production (Chen2006) - obtained from fig 6 & equ 17 & 18
    parser.addParameter('O2_n', 200);           % [uM] ; tissue O2 
    % concentration in the neuron (M.E.)
    parser.addParameter('K_mO2_n', 243);        % [uM] ; Chen2006
    parser.addParameter('LArg_n', 100);         % [uM] ; 
    parser.addParameter('K_mArg_n', 1.5);       % [uM] ;     
    parser.addParameter('k_O2_n', 9.6e-9);      % [uM^-2 ms^-1] ; % (Kavdia2002)
    parser.addParameter('x_nk', 25);            % [um] ;  (M.E.)
    parser.addParameter('D_cNO', 3300e-3);      % [um^2 ms^-1] ; Diffusion 
    % coefficient NO (Malinski1993)    
    parser.addParameter('V_spine', 8e-5);       % [nL] ; volume of the 
    % neuronal dendritic spine Santucci2008
    parser.addParameter('k_ex', 1600e-3);       % [ms^-1] ; decay rate constant
    % of internal calcium concentration Santucci2008
    parser.addParameter('Ca_rest', 0.1);        % [uM] ; resting calcium 
    % concentration (in Comerford+David2008: 2.830 mM; in Santucci2008P: 0.1 \muM)
    parser.addParameter('lambda_buf', 20);      % [-] ; buffer capacity Santucci2008
    parser.addParameter('V_maxNOS', 25e-6);     % [uM ms^-1] ; M.E.
    parser.addParameter('K_actNOS', 9.27e-2);   % [uM] ; 
    parser.addParameter('mu2_n', 0.0167e-3);    % [ms^-1] ; rate constant at
    % which the nNOS is deactivated Comerford2008

    % input 
    parser.addParameter('startpulse', 100);    
    parser.addParameter('lengthpulse', 16);
    parser.addParameter('Whiskerpad', 1);
    parser.addParameter('secondpulse', 124);
    parser.addParameter('secondlength', 1);
    
    parser.addParameter('dt', 0.001);
    parser.addParameter('XLIM2', 150);
    parser.addParameter('input_data', linspace(0, 1200, 1000)); 
    
    parser.parse(varargin{:})
    params = parser.Results;
end

function u0 = initial_conditions(idx)
    u0 = zeros(length(fieldnames(idx)), 1);

    u0(idx.E_t) = 0.0888;       % [-]
    u0(idx.I_t) = 0.3330;       % [-]
    u0(idx.K_e) = 3.5;          % [mM]
    u0(idx.Na_sa) = 9.37;       % [mM]
    u0(idx.Na_d) = 9.42;        % [mM]
    u0(idx.O2) = 0.0266;        % [uM]
    u0(idx.CBV) = 1.3167;       % [-]
    u0(idx.HbR) = 0.6665;       % [-]
    u0(idx.Ca_n) = 0.1;         % [uM]
    u0(idx.nNOS_act_n) = 0.318; % [uM]
    u0(idx.NO_n) = 0.1671;      % [uM]
end