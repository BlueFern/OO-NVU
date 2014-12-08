classdef NVU < handle
    properties
        astrocyte
        wall
        smcec
        i_astrocyte
        i_wall
        i_smcec
        n
    end
    methods 
        function self = NVU(varargin)
            self.astrocyte = Astrocyte();
            self.wall = WallMechanics();
            self.smcec = SMCEC('J_PLC', 0.4);
            na = length(fieldnames(self.astrocyte.index));
            nw = length(fieldnames(self.wall.index));
            ns = length(fieldnames(self.smcec.index));
            self.i_astrocyte = 1:na;
            self.i_smcec = na + (1:ns);
            self.i_wall = na + ns + (1:nw);
            self.n = na + ns + nw;
        end
        function du = rhs(self, t, u)
            % Separate out the model components
            ua = u(self.i_astrocyte, :);
            us = u(self.i_smcec, :);
            uw = u(self.i_wall, :);
            
            % Evaluate the coupling pieces needed for RHS evaluation
            K_p = self.astrocyte.output(t, ua);
            [J_KIR_i, Ca_i] = self.smcec.output(t, us, K_p);
            [R, h] = self.wall.output(t, uw);
            
            du = zeros(size(u));
            du(self.i_astrocyte, :) = self.astrocyte.rhs(t, ua, J_KIR_i);
            du(self.i_wall, :) = self.wall.rhs(t, uw, Ca_i);
            du(self.i_smcec, :) = self.smcec.rhs(t, us, R, h, K_p);
        end
        function out = u0(self)
            out = zeros(self.n, 1);
            out(self.i_astrocyte) = self.astrocyte.u0;
            out(self.i_smcec) = self.smcec.u0;
            out(self.i_wall) = self.wall.u0;
        end
        
    end
end