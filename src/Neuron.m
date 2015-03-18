classdef Neuron < handle
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
        function [du, varargout] = rhs(self, t, u)
            % Initalise inputs and parameters
            t = t(:).';
            p = self.params;
            idx = self.index;
            
            % State variables
            N_Na_n = u(idx.N_Na_n, :);
            N_K_n = u(idx.N_K_n, :);
            
            Ca_n = u(idx.Ca_n, :);                  % NO pathway
            nNOS_act_n = u(idx.nNOS_act_n, :);
            NO_n = u(idx.NO_n, :);
            
            % Algebraic equations            
            J_Na_n = p.k_C * self.input_f(t);
            J_K_n = -J_Na_n;
            
            w_NR2A = self.input_Glu(t) ./ (p.K_m_A + self.input_Glu(t)); 
            w_NR2B = self.input_Glu(t) ./ (p.K_m_B + self.input_Glu(t)); 
            I_Ca = (-4*p.v_n*p.G_M*p.P_Ca_P_M*(p.Ca_ex/p.M))/(1+exp(-0.08*(v_n+20)))...
                   *(exp(2*1e-3*v_n*Farad/(R*T)))/(1-exp(2*1e-3*v_n*Farad/(R*T)));
            I_Ca_tot = I_Ca *(0.63*NE(flu.P_NR2AO)+11*NE(flu.P_NR2BO));
            phi_N = 
            dphi_N =
            m_c =
            CaM = 
            
            
            % Differential equations:
            du(idx.N_Na_n, :) = J_Na_n;
            du(idx.N_K_n, :) = -du(idx.N_Na_n, :);
            du(idx.Ca_n
            du(idx.nNOS_act
            du(idx.NO_n
            du(idx.
            du(idx.
            
            du = bsxfun(@times, self.enabled, du);
            if nargout == 2
               Uout = zeros(self.n_out, size(u, 2));
               Uout(self.idx_out.ft, :) = self.input_f(t);
               Uout(self.idx_out.J_Na_n, :) = J_Na_n;
               Uout(self.idx_out.J_K_n, :) = J_K_n;
               varargout = {Uout};
            end
        end
        function f = input_f(self, t)      
            p = self.params;                
            f = zeros(size(t));
            ii = p.t_0 <= t & t < p.t_1;
            f(ii) = ...
                p.F_input * p.gab / ...
                (p.ga * p.gb) * ...
                (1 - (t(ii) - p.t_0) / p.delta_t).^(p.beta - 1) .* ...
                ((t(ii) - p.t_0) / p.delta_t).^(p.alpha - 1);
            f(p.t_2 <= t & t <= p.t_3) = -p.F_input;
        end
        function J_Na_n = shared(self, t)    %shared variable
            t = t(:).';
            p = self.params;
            J_Na_n = p.k_C * self.input_f(t);
        end
        function Glu = input_Glu(self, t)
            p = self.params;
            Glu = (p.Glu_max - p.Glu_min) * ( ...
            0.5 * tanh((t - p.t_0_Glu) / p.theta_L_Glu) - ...
            0.5 * tanh((t - p.t_2_Glu) / p.theta_R_Glu)) + p.Glu_min;
        end
        function names = varnames(self)
            names = [fieldnames(self.index); fieldnames(self.idx_out)];
        end
    end
end    
        
function idx = indices()    %for variables with odes
idx.N_Na_n = 1;
idx.N_K_n = 2;
end        

function [idx, n] = output_indices()    %for variables in nargout loop
idx.ft = 1;
idx.J_Na_n = 2;
idx.J_K_n = 3; 
n = numel(fieldnames(idx));
end
        
function params = parse_inputs(varargin)
parser = inputParser();

% input 
parser.addParameter('startpulse', 200);     %Used for f(t)
parser.addParameter('lengthpulse', 200);
parser.addParameter('lengtht1', 10);
parser.addParameter('F_input', 2.5); %s  
parser.addParameter('alpha', 2);
parser.addParameter('beta', 5);
parser.addParameter('delta_t', 10); %s
parser.addParameter('k_C', 7.35e-5); %uM m s^-1
parser.addParameter('startpulse_Glu', 200); % Glu input 
parser.addParameter('lengthpulse_Glu', 200); % Glu input 
parser.addParameter('Glu_max', 1846); % microM (one vesicle, Santucci2008)
parser.addParameter('Glu_min', 0); % microM
parser.addParameter('theta_L_Glu', 1); % slope of Glu input 
parser.addParameter('theta_R_Glu', 1); % slope of Glu input 

parser.addParameter('K_m_A', 650); % [uM] - fit to Santucci2008
parser.addParameter('K_m_B', 2800); % [uM] - fit to Santucci2008

parser.addParameter('v_n', -40);          % [mV] ; the neuronal membrane potential , assumed to be approx constant in this model
parser.addParameter('G_M', 46);        	% [pS] ; the conductance of the NMDA channel to Ca2+ compaired  
parser.addParameter('P_Ca_P_M', 3.6);         	% [-] ; the relative conductance of the NMDA channel to Ca2+ compared to monovalent ions
parser.addParameter('Ca_ex', 2e3);          % [microM] ; the external calcium concentration (in Comerford+David2008: 1.5 mM!)
parser.addParameter('M', 1.3e5);        % [microM] ; the concentration of monovalent ions in the neuron
parser.addParameter('R_gas', 8.314);      	% [Jmol^{-1}K^{-1}]
parser.addParameter('T', 310.65);     	% [K] ; temperature


parser.addParameter('Q1', 1.9e5);      	% [-]
parser.addParameter('Q2', 2.1e5);      	% [-]
parser.addParameter('Q3', 0.4e5);      	% [-]
parser.addParameter('Q4', 0.26e5);      	% [-]

parser.parse(varargin{:})
params = parser.Results;
params.t_0 = params.startpulse;
params.t_1 = params.t_0 + params.lengtht1;
params.t_2 = params.t_0 + params.lengthpulse;
params.t_3 = params.t_1 + params.lengthpulse;
params.gab = factorial(params.alpha + params.beta - 1);
params.ga = factorial(params.alpha - 1);
params.gb = factorial(params.beta - 1);
params.t_0_Glu = params.startpulse_Glu;
params.t_2 = params.t_0_Glu + params.lengthpulse_Glu;
end

function u0 = initial_conditions(idx)
u0 = zeros(length(fieldnames(idx)), 1);

%ICs for neuron? Set so that initial Na is 0 and K goes to 0. Need
%physiological values
u0(idx.N_Na_n) = 0;
u0(idx.N_K_n) = 1.8372353094e-3;
end






