classdef NVU_red < handle
    % The 'NVU' Code wraps the three subcodes (Astrocyte, SMCEC and Wall
    %   Mechanics) together to be solved by the ode15s solver for stiff
    %   problems. 
    properties
        wall
        smcec
        i_wall
        i_smcec
        offsets
        outputs
        n
        u0
        T
        U
        odeopts
        elapsedtime
    end
    methods 
        function self = NVU_red(wall, smcec, varargin)
            % Deal with input parameters
            params = parse_inputs(varargin{:});
            names = fieldnames(params);
            for i = 1:numel(names)
                self.(names{i}) = params.(names{i});
            end
            self.wall = wall;
            self.smcec = smcec;
            
            % Construct mapping to full state vector
            nw = length(fieldnames(self.wall.index));
            ns = length(fieldnames(self.smcec.index));
            self.offsets = [0, nw];
            self.outputs = {[],[]};
            self.i_wall = 1:nw;
            self.i_smcec = nw + (1:ns);
            self.n = nw + ns;
            
            
            
            self.init_conds()
        end
        function du = rhs(self, t, u)
            % Separate out the model components
            us = u(self.i_smcec, :);
            uw = u(self.i_wall, :);
            
            % Evaluate the coupling quantities to be passed between
            % submodels as coupling
            [Ca_i] = self.smcec.shared(t, us);
            [R, h] = self.wall.shared(t, uw);
            
            du = zeros(size(u));
            du(self.i_wall, :) = self.wall.rhs(t, uw, Ca_i);
            du(self.i_smcec, :) = self.smcec.rhs(t, us, R, h);
        end
        function init_conds(self)
            self.u0 = zeros(self.n, 1);
            self.u0(self.i_smcec) = self.smcec.u0;
            self.u0(self.i_wall) = self.wall.u0;
        end
        function simulate(self)
            self.init_conds()
            f = @(t, u) self.rhs(t, u);
            tic
            [self.T, self.U] = ode15s(f, self.T, self.u0, self.odeopts);
            % Now evaluate all of the additional parameters
            us = self.U(:, self.i_smcec).';
            uw = self.U(:, self.i_wall).';
            
            [Ca_i] = self.smcec.shared(self.T, us);
            [R, h] = self.wall.shared(self.T, uw);
            
            [~, self.outputs{1}] = self.wall.rhs(self.T, uw, Ca_i);
            [~, self.outputs{2}] = self.smcec.rhs(self.T, us, R, h);

            
            self.elapsedtime = toc;
        end
        function u = out(self, input_str)
            success = false;
            modules = {self.wall, self.smcec};
            for i = 1:2
                module = modules{i};
                if ismember(input_str, fieldnames(module.index))
                   u = self.U(:, self.offsets(i) + (module.index.(input_str)));
                   success = true;
                   break
                end
                if ismember(input_str, fieldnames(module.idx_out))
                    u = self.outputs{i}(module.idx_out.(input_str), :);
                    success = true;
                    break
                end
            end
            if ~success
                error('NVU:InvalidFieldName', 'No matching field: %s', input_str)
            end
        end
        
    end
end

function params = parse_inputs(varargin)
parser = inputParser();
parser.addParameter('odeopts', odeset());
parser.addParameter('T', linspace(0, 600, 1000));
parser.parse(varargin{:});
params = parser.Results;
end