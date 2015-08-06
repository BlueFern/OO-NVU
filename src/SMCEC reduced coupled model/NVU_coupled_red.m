classdef NVU_coupled_red < handle
    % The 'NVU_coupled' Code wraps the 6 subcodes (Astrocyte, SMCEC and Wall
    %   Mechanics for both systems) together to be solved by the ode15s 
    %   solver for stiff problems. 
    
    
    properties
        wall_1
        wall_2
        smcec_1
        smcec_2
        i_wall_1
        i_wall_2
        i_smcec_1
        i_smcec_2
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
        function self = NVU_coupled_red(wall_1, wall_2, smcec_1, smcec_2, varargin)
            % Deal with input parameters
            params = parse_inputs(varargin{:});
            names = fieldnames(params);
            for i = 1:numel(names)
                self.(names{i}) = params.(names{i});
            end
            self.wall_1 = wall_1;
            self.wall_2 = wall_2;
            self.smcec_1 = smcec_1;
            self.smcec_2 = smcec_2;
            
            % Construct mapping to full state vector
            nw_1 = length(fieldnames(self.wall_1.index));
            nw_2 = length(fieldnames(self.wall_2.index));
            ns_1 = length(fieldnames(self.smcec_1.index));
            ns_2 = length(fieldnames(self.smcec_2.index));
            self.offsets = [0, nw_1, nw_1+nw_2, nw_1+nw_2+ns_1, nw_1+nw_2+ns_1+ns_2];
            self.outputs = {[],[],[],[]};
            self.i_wall_1 = 1:nw_1;
            self.i_wall_2 = nw_1 + (1:nw_2);
            self.i_smcec_1 = nw_1 + nw_2 + (1:ns_1);
            self.i_smcec_2 = nw_1 + nw_2 + ns_1 + (1:ns_2);
            self.n = nw_1 + nw_2 + ns_1 + ns_2;

            self.init_conds()
        end
        function du = rhs(self, t, u)
            % Separate out the model components
            us_1 = u(self.i_smcec_1, :);
            us_2 = u(self.i_smcec_2, :);
            uw_1 = u(self.i_wall_1, :);
            uw_2 = u(self.i_wall_2, :);
            
            % Evaluate the coupling quantities to be passed between
            % submodels as coupling
            [Ca_i_1, I_i_1, v_i_1, Ca_j_1, I_j_1, v_j_1] = self.smcec_1.shared(t, us_1);
            [Ca_i_2, I_i_2, v_i_2, Ca_j_2, I_j_2, v_j_2] = self.smcec_2.shared(t, us_2);
            [R_1, h_1] = self.wall_1.shared(t, uw_1);
            [R_2, h_2] = self.wall_2.shared(t, uw_2);
            
            du = zeros(size(u));
            du(self.i_wall_1, :) = self.wall_1.rhs(t, uw_1, Ca_i_1);
            du(self.i_wall_2, :) = self.wall_2.rhs(t, uw_2, Ca_i_2);
            du(self.i_smcec_1, :) = self.smcec_1.rhs(t, us_1, R_1, h_1, Ca_i_2, I_i_2, v_i_2, Ca_j_2, I_j_2, v_j_2);
            du(self.i_smcec_2, :) = self.smcec_2.rhs(t, us_2, R_2, h_2, Ca_i_1, I_i_1, v_i_1, Ca_j_1, I_j_1, v_j_1);
        end
        function init_conds(self)
            self.u0 = zeros(self.n, 1);
            self.u0(self.i_smcec_1) = self.smcec_1.u0;
            self.u0(self.i_smcec_2) = self.smcec_2.u0;
            self.u0(self.i_wall_1) = self.wall_1.u0;
            self.u0(self.i_wall_2) = self.wall_2.u0;
        end
        function simulate(self)
            self.init_conds()
            f = @(t, u) self.rhs(t, u);
            tic
            [self.T, self.U] = ode15s(f, self.T, self.u0, self.odeopts);
            % Now evaluate all of the additional parameters
            us_1 = self.U(:, self.i_smcec_1).';
            us_2 = self.U(:, self.i_smcec_2).';
            uw_1 = self.U(:, self.i_wall_1).';
            uw_2 = self.U(:, self.i_wall_2).';
            
            [Ca_i_1, I_i_1, v_i_1, Ca_j_1, I_j_1, v_j_1] = self.smcec_1.shared(self.T, us_1);
            [Ca_i_2, I_i_2, v_i_2, Ca_j_2, I_j_2, v_j_2] = self.smcec_2.shared(self.T, us_2);
            [R_1, h_1] = self.wall_1.shared(self.T, uw_1);
            [R_2, h_2] = self.wall_2.shared(self.T, uw_2);

            [~, self.outputs{1}] = self.wall_1.rhs(self.T, uw_1, Ca_i_1);
            [~, self.outputs{2}] = self.wall_2.rhs(self.T, uw_2, Ca_i_2);
            [~, self.outputs{3}] = self.smcec_1.rhs(self.T, us_1, R_1, h_1, Ca_i_2, I_i_2, v_i_2, Ca_j_2, I_j_2, v_j_2);
            [~, self.outputs{4}] = self.smcec_2.rhs(self.T, us_2, R_2, h_2, Ca_i_1, I_i_1, v_i_1, Ca_j_1, I_j_1, v_j_1);
            
            self.elapsedtime = toc;
        end
        function u = out(self, input_str)
            success = false;
            modules = {self.wall_1, self.wall_2, self.smcec_1, self.smcec_2};
            for i = 1:4
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