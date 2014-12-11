%% WallMechanics MATLAB model documentation
% This page describes the model components present in the Wall Mechanics
% model.

%% Right hand side function
% The RHS function for the |WallMechanics| model is of the following form
%%
% 
%   function du = rhs(self, t, u, Ca_i)
% 
%%
% This form shows which model components are required to couple with the
% |WallMechanics| model. These components are generated by calling the
% |shared| method of the other model components. Currently these are |Ca_i|
% from the |SMCEC| model.


%% Sharing variables
% The |shared| method generates the coupling variables that need to be sent
% to other components of the NVU model. The function form is
%%
%
%     function [R, h] = shared(self, ~, u)
%
%%
% If additional variables are required to pass to other model components,
% this is the method that you will need to edit.

A = WallMechanics();

%% State variables
% The following state variables are present in the |WallMechanics| model. These
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
% The |WallMechanics| model has the following parameters, with specified default
% values. To change these, these can be called at object creation with, for
% example
%
%     W = WallMechanics('eta', 2e5)
%
% To change them later, they can simply be modified by, for example:
%    
%     W.params.eta = 2e5;
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