classdef WallMechanics_1 < handle
    % The 'WallMechanics' code contains the following sections of the model:
    %   The Contraction and Mechanical Models.
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
        function self = WallMechanics_1(varargin)
            self.params = parse_inputs(varargin{:});
            self.index = indices();
            self.u0 = initial_conditions(self.index);
            self.enabled = true(size(self.u0));
            [self.idx_out, self.n_out] = output_indices();
        end
        function [du, varargout] = rhs(self, t, u, Ca_i_1)
            % Initalise inputs and parameters
            idx = self.index;
            p = self.params;
            % Make the output function compute these so as to avoid
            % duplicating equations
            [R_1, h_1] = self.shared(t, u);
            
            Mp_1 = u(idx.Mp_1, :);
            AMp_1 = u(idx.AMp_1, :);
            AM_1 = u(idx.AM_1, :);
            
            %% Contraction Equations
            K_1_1 = p.gamma_cross * Ca_i_1.^3;
            K_6_1 = K_1_1;
            M_1 = 1 - AM_1 - AMp_1 - Mp_1;
            du(idx.Mp_1, :) = p.K_4 * AMp_1 + K_1_1 .* M_1 - (p.K_2 + p.K_3) * Mp_1;
            du(idx.AMp_1, :) = p.K_3 * Mp_1 + K_6_1 .* AM_1 - (p.K_4 + p.K_5) * AMp_1;
            du(idx.AM_1, :) = p.K_5 * AMp_1 - (p.K_7 + K_6_1) .* AM_1;
            
            % Mechanical Equations
            F_r_1 = AMp_1 + AM_1;
            E_1 = p.E_passive + F_r_1 * (p.E_active - p.E_passive);
            R_0_1 = p.R_0_passive + F_r_1 * (p.alpha - 1) * p.R_0_passive;
            du(idx.R_1, :) = p.R_0_passive / p.eta * ( R_1 * p.P_T ./ h_1 - ...
                E_1 .* (R_1 - R_0_1) ./ R_0_1);
            
            du = bsxfun(@times, self.enabled, du);
            
            if nargout == 2
               Uout = zeros(self.n_out, size(u, 2));
               Uout(self.idx_out.M_1, :) = M_1;
               Uout(self.idx_out.F_r_1, :) = F_r_1;
               varargout{1} = Uout;
            end
        end
        function [R_1, h_1] = shared(self, ~, u)
            R_1 = u(self.index.R_1, :);
            h_1 = 0.1 * R_1;
        end     
        function names = varnames(self)
            names = [fieldnames(self.index); fieldnames(self.idx_out)];
        end
    end 
end

function idx = indices()
% Index of parameters needing inital conditions 
idx.Mp_1 = 1;
idx.AMp_1 = 2;
idx.AM_1 = 3;
idx.R_1 = 4;
end

function [idx, n] = output_indices()
% Index of all other output parameters
idx.M_1 = 1;
idx.F_r_1 = 2;
n = numel(fieldnames(idx));
end

function params = parse_inputs(varargin)
parser = inputParser();
% Contraction Equation Constants
parser.addParameter('K_2', 0.5); % s^-1
parser.addParameter('K_3', 0.4); % s^-1
parser.addParameter('K_4', 0.1); % s^-1
parser.addParameter('K_5', 0.5); % s^-1
parser.addParameter('K_7', 0.1); % s^-1
parser.addParameter('gamma_cross', 17); %uM^-3 s^-1
% Mechanical Equation Constants
parser.addParameter('eta', 1e4); %Pa s
parser.addParameter('R_0_passive', 20e-6); % m
parser.addParameter('h_0_passive', 3e-6); % m
parser.addParameter('P_T', 4000); % Pa
parser.addParameter('E_passive', 66e3); % Pa
parser.addParameter('E_active', 233e3); % Pa
parser.addParameter('alpha', 0.6); % [-]

parser.parse(varargin{:});
params = parser.Results;

end
function u0 = initial_conditions(idx)
u0 = zeros(length(fieldnames(idx)), 1);
% Inital estimations of parameters from experimental data
u0(idx.Mp_1) = .25;
u0(idx.AMp_1) = .25;
u0(idx.AM_1) = .25;
u0(idx.R_1) = 15e-6;
end
