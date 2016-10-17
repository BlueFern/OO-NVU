%% Demonstration script for new NVU model version 1.2

%% Version control
% This is the latest version of the NVU model. It includes:
% TRPV4
% NO pathway
% Astrocytic Ca2+

%% Construct NVU
% The NVU consists of a number of submodules, implemented as MATLAB
% classes, presently an astrocyte, a lumped SMC/EC model, and a model of
% the wall mechanics. 
%
% The parameters of each of these submodules are
% specified when the modules are constructed, here in the call to NVU, see
% for example the |SMCEC| part of the |NVU| call below:
%
% Options for the ODE solver (currently |ode15s|) are provided by
% specifying the |odeopts| parameter. The code works fine with default
% tolerances.
clear all
clc
odeopts = odeset('RelTol', 1e-04, 'AbsTol', 1e-04, 'MaxStep', 0.5, 'Vectorized', 1);

% Things to change:
XLIM1 = 0; XLIM2 = 400;   % x limits of any plots
PRESSURE = 4000;          % Pressure in Pa, change here!! Because the parameter is used in both SMCEC and WallMechanics ***********
START_PULSE  = 100;      % Start time of neuronal stimulation
LENGTH_PULSE = 100;      % Length of stimulation 
GLU_SWITCH   = 1;        % Turn on glutamate input to SC, default 1
NO_SWITCH    = 1;        % Turn on Nitric Oxide production, default 1
TRPV4_SWITCH = 0;        % Turn on TRPV4 channel flux in AC, default 1
J_PLC        = 0.18;     % Jplc value in EC: 0.18 for steady state (default), 0.4 for oscillations (vasomotion).

% Initiate the NVU
nv = NVU(Neuron('startpulse', START_PULSE, 'lengthpulse', LENGTH_PULSE, 'GluSwitch', GLU_SWITCH, 'NOswitch', NO_SWITCH), ...
    Astrocyte('startpulse', START_PULSE, 'lengthpulse', LENGTH_PULSE, 'rhoSwitch', GLU_SWITCH, 'TRPV4switch', TRPV4_SWITCH), ...
    WallMechanics('P_T', PRESSURE), ...
    SMCEC('J_PLC', J_PLC, 'NOswitch', NO_SWITCH, 'P_T', PRESSURE), 'odeopts', odeopts);

nv.T = linspace(0, XLIM2, 5000);     % Initiate time vector
nv.simulate();                       % Run simulation

%% Plot stuff
% Either run plot_all_variables.m to get everything, or write your own code
% here to plot specific variables or fluxes. 
% E.g.

figure(10);         % Initialise figure

subplot(1,2,1)
plot(nv.T, nv.out('J_stretch_i'));
xlabel('Time [s]'); title('Ca2+ flux through SMC stretch channel [\muM/s]')
xlim([XLIM1 XLIM2])

subplot(1,2,2)
plot(nv.T, nv.out('J_stretch_j'))
xlabel('Time [s]'); title('Ca2+ flux through EC stretch channel [\muM/s]')
xlim([XLIM1 XLIM2])

%%
figure(11);
plot(nv.T, nv.out('s_k'));

