classdef WallMechanics < handle
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
        function self = WallMechanics(varargin)
            self.params = parse_inputs(varargin{:});
            self.index = indices();
            self.u0 = initial_conditions(self.index);
            self.enabled = true(size(self.u0));
            [self.idx_out, self.n_out] = output_indices();
        end
        function [du, varargout] = rhs(self, t, u, Ca_i, R_cGMP2)
            % Initalise inputs and parameters
            idx = self.index;
            p = self.params;
            % Make the output function compute these so as to avoid
            % duplicating equations
            [R, h] = self.shared(t, u);
            [pR] = self.inflate(t);
            Mp = u(idx.Mp, :);
            AMp = u(idx.AMp, :);
            AM = u(idx.AM, :);
            
            %% Contraction Equations
            K_1 = p.gamma_cross * Ca_i.^p.n_cross;
            K_6 = K_1;
            
%             K_2 = 0.5; % NO excluded
            K_2 = 58.1395 * p.k_mlcp_b + 58.1395 * p.k_mlcp_c * R_cGMP2;  % 17.64 / 16.75 errechnet sich aus dem Shift, um K2_c = 0.5 bei der baseline zu bekommen, NO included - muss vllt noch geaendert werden! 
            K_5 = K_2;
            
            M = 1 - AM - AMp - Mp;
            du(idx.Mp, :) = p.K_4 * AMp + K_1 .* M - (K_2 + p.K_3) .* Mp;
            du(idx.AMp, :) = p.K_3 * Mp + K_6 .* AM - (p.K_4 + K_5) .* AMp;
            du(idx.AM, :) = K_5 .* AMp - (p.K_7 + K_6) .* AM;
            
            % Mechanical Equations
            F_r = AMp + AM;
            E = p.E_passive + F_r * (p.E_active - p.E_passive);
            R_0 = p.R_0_passive + F_r * (p.alpha - 1) * p.R_0_passive;
           % R_new = self.input_R(t);

            du(idx.R, :) =  p.R_0_passive / p.eta * ( R * pR .* p.P_T ./ h - ...
                            E .* (R - R_0) ./ R_0);
            
            du = bsxfun(@times, self.enabled, du);
            
            if nargout == 2
               Uout = zeros(self.n_out, size(u, 2));
               Uout(self.idx_out.F_r, :) = F_r;
               Uout(self.idx_out.K_2, :) = K_2;
               Uout(self.idx_out.P_T_t_wall, :) = p.P_T;
               Uout(self.idx_out.M, :) = M;
               varargout{1} = Uout;
            end
        end
%         function R_input = input_R(~, t)
%             % Input signal; the smooth pulse function rho
%             R_input = 7.74e-6* (0.5* tanh((t-110)/0.8)-0.5* tanh((t-115)/1.9))+19.37e-6; (Joerik)
%         end
        function [pR] = inflate(self, t)
%             t_scalar = t(end,:);
            pR = 1; %5*(0.5 .* tanh((t-75)./0.8)-0.5.* tanh((t-76) ./1.9)); %u(self.index.R, :); (Joerik)
            
        end   
        function [R, h] = shared(self, ~,u)
           
            R = u(self.index.R, :);
            h = 0.1 * R;
            
        end  
        function names = varnames(self)
            names = [fieldnames(self.index); fieldnames(self.idx_out)];
        end
        
        function Pressurechange = input_Pressure(self, t) %TimvdBoom
            
            p = self.params;
            Pressurechange = 0.0075*p.PressureSwitch * (p.Pressure_change - p.P_T) * ( ...
                0.5 * tanh((t - p.Pressurechange_t_1) / 10) - ...
                0.5 * tanh((t - p.Pressurechange_t_2) / 10)) + 0.0075*p.P_T;
        end
        
    end 
end

function idx = indices()
    % Index of parameters needing inital conditions 
    idx.Mp = 1;
    idx.AMp = 2;
    idx.AM = 3;
    idx.R = 4;
end

function [idx, n] = output_indices()
    % Index of all other output parameters
    idx.M = 1;
    idx.F_r = 2;
    idx.R = 3;
    idx.K_2 = 4;
    idx.P_T_t_wall=5;
    n = numel(fieldnames(idx));
end

function params = parse_inputs(varargin)
    parser = inputParser();
    % Contraction Equation Constants
    % parser.addParameter('K_2', 0.5); % s^-1
    parser.addParameter('K_3', 0.4); % s^-1
    parser.addParameter('K_4', 0.1); % s^-1
    % parser.addParameter('K_5', 0.5); % s^-1
    parser.addParameter('K_7', 0.1); % s^-1
    parser.addParameter('gamma_cross', 17); %uM^-3 s^-1
    parser.addParameter('n_cross', 3); % fraction constant of the phosphorylation crossbridge
    % Mechanical Equation Constants
    parser.addParameter('eta', 1e4); %Pa s
    parser.addParameter('R_0_passive', 20e-6); % m
    parser.addParameter('P_T', 4000); % Pa
    parser.addParameter('E_passive', 66e3); % Pa
    parser.addParameter('E_active', 233e3); % Pa
    parser.addParameter('alpha', 0.6); % [-]
    parser.addParameter('k_mlcp_b', 0.0086); % [s^-1]
    parser.addParameter('k_mlcp_c', 0.0327); % [s^-1]
    
    parser.addParameter('PressureSwitch',1);
    parser.addParameter('Pressure_change', 4000); % Pressure in Pa
    
    parser.parse(varargin{:});
    params = parser.Results;
    params.Pressurechange_t_1=1000;
    params.Pressurechange_t_2=2000;

end
function u0 = initial_conditions(idx)
    u0 = zeros(length(fieldnames(idx)), 1);
    % Inital estimations of parameters from experimental data
    u0(idx.Mp) = .25;
    u0(idx.AMp) = .25;
    u0(idx.AM) = .25;
    u0(idx.R) = 24.8e-6; %15e-6;
end
