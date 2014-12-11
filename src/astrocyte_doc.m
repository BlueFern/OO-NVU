%% Astrocyte model documentation
% This page describes the model components present in the astrocyte model. 

%% Right hand side function
% The RHS function for the astrocyte is of the following form
%%
% 
%   function du = rhs(self, t, u, J_KIR_i)
% 
%%
% The important thing is that the astrocyte model requires J_KIR_i as a
% coupling variable. Note that currently the neuron model is incorporated
% within this model, it should in time be broken out as a separate class.

%% Sharing variables
% The |shared| method generates the coupling variables that need to be sent
% to other components of the NVU model. The function form is
%%
%
%     function K_p = shared(self, ~, u)
%
%%
% If additional variables are required to pass to other model components,
% this is the method that you will need to edit.

A = Astrocyte();

%% State variables
% The following state variables are present in the Astrocyte model. These
% can be retrieved from the NVU after simulation by doing |nv.out(varname)|
% where |nv| is an NVU object
vars = fieldnames(A.index);
for i = 1:numel(vars)
    fprintf('%s\n', vars{i});
end

%% Additional quantities
% The following additional quantities can also be retrieved in the same way
% as the state variables after simulation. To add more of these, add more
% entries in the |output_indices| method, together with an associated entry
% at the end of the |rhs| method (inside the |if nargout == 2| clause).
vars = fieldnames(A.idx_out);
for i = 1:numel(vars)
    fprintf('%s\n', vars{i});
end

%% Parameters
% The astrocyte model has the following parameters, with specified default
% values. To change these, these can be called at object creation with, for
% example
%
%     A = Astrocyte('g_K_k', 100)
%
% To change them later, they can simply be modified by, for example:
%    
%     A.params.g_K_k = 100;
paramnames = fieldnames(A.params);
for i = 1:numel(paramnames)
   fprintf('%-15s%g\n', strcat(paramnames{i}, ' :'), A.params.(paramnames{i}))  
end

%% Initial conditions
% The modeled initial conditions are as follows:
vars = fieldnames(A.index);
for i = 1:numel(vars)
    fprintf('%-15s%g\n', strcat(vars{i}, ' :'), A.u0(A.index.(vars{i})));
end




