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
            t = t(:).';
            p = self.params;
            idx = self.index;
            
            N_Na_n = u(idx.N_Na_n, :);
            N_K_n = u(idx.N_K_n, :);
            
            J_Na_n = p.k_C * self.input_f(t);
            J_K_n = -J_Na_n;
            
            du(idx.N_Na_n, :) = J_Na_n;
            du(idx.N_K_n, :) = -du(idx.N_Na_n, :);

            
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
        function J_Na_n = shared(self, t, f)    %shared variable
            t = t(:).';
            p = self.params;
            J_Na_n = p.k_C * self.input_f(t);
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

parser.addParameter('startpulse', 200);     %Used for f(t)
parser.addParameter('lengthpulse', 200);
parser.addParameter('lengtht1', 10);
parser.addParameter('F_input', 2.5); %s  
parser.addParameter('alpha', 2);
parser.addParameter('beta', 5);
parser.addParameter('delta_t', 10); %s
parser.addParameter('k_C', 7.35e-5); %uM m s^-1

parser.parse(varargin{:})
params = parser.Results;
params.t_0 = params.startpulse;
params.t_1 = params.t_0 + params.lengtht1;
params.t_2 = params.t_0 + params.lengthpulse;
params.t_3 = params.t_1 + params.lengthpulse;
params.gab = factorial(params.alpha + params.beta - 1);
params.ga = factorial(params.alpha - 1);
params.gb = factorial(params.beta - 1);
end

function u0 = initial_conditions(idx)
u0 = zeros(length(fieldnames(idx)), 1);

%ICs for neuron? Set so that initial Na is 0 and K goes to 0. Need
%physiological values
u0(idx.N_Na_n) = 0;
u0(idx.N_K_n) = 1.8372353094e-3;
end






