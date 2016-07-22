%% Demonstration script for new NVU model

% Version control
% This is the latest version of the NVU model. It includes:
% - basic K+ / Ca2+ 

% Construct NVU
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

% Handy parameters to modify
START_PULSE  = 200;      % Time of neuronal stimulation
LENGTH_PULSE = 40;       % Length of stimulation (used if CURREN_TYPE = 2)
GLU_SWITCH   = 1;        % Turn on glutamate input to SC
NAT_SWITCH   = 1;        % Turn on NaT channel in dendrite
BUFF_SWITCH  = 1;        % Turn on K+ buffering of ECS by AC
NO_SWITCH    = 1;        % Turn on Nitric Oxide production 
CURRENT_TYPE = 1;        % 1 for Gaussian, 2 for tanh function
J_PLC        = 0.18;     % Jplc value in EC: 0.18 for steady state, 0.4 for oscillations

nv = NVU(Neuron('startpulse', START_PULSE, 'lengthpulse', LENGTH_PULSE, 'GluSwitch', GLU_SWITCH, 'NaTswitch', NAT_SWITCH, 'buffSwitch', BUFF_SWITCH, 'NOswitch', NO_SWITCH, 'currentType',CURRENT_TYPE), ...
    Astrocyte('startpulse', START_PULSE, 'lengthpulse', LENGTH_PULSE, 'rhoSwitch', GLU_SWITCH), ...
    WallMechanics(), ...
    SMCEC('J_PLC', J_PLC, 'NOswitch', NO_SWITCH), 'odeopts', odeopts);

nv.T = linspace(0, 400, 10000);    
    
nv.simulate()

% legend('no NaT, no buffer','NaT, no buffer','no NaT, buffer','NaT,buffer')

% Lists of quantities that can be retrieved from the NVU model after
% simulation are given in the documentation pages of the individual model
% components.
%
% Model components are retrieved using the |out| method as follows, where
% we do a simple plot of intracellular calcium from the SMC:

% Plot, e.g. Ca_i
% plot(nv.T, nv.out('Ca_i'))
% xlabel('time (s)')
% ylabel('[Ca^{2+}] (\muM)')

figure(6);
subplot(3,2,1)
hold all;
plot(nv.T, nv.out('R'))
xlabel('Time [s]'); ylabel('R')
xlim([100 400])

subplot(3,2,2)
hold all;
plot(nv.T, nv.out('K_s')/1e3)
xlabel('Time [s]'); ylabel('K_s')
xlim([100 400])

subplot(3,2,3)
hold all;
plot(nv.T, nv.out('v_sa'));
xlabel('Time [s]'); ylabel('v_{sa}')
xlim([100 400])

subplot(3,2,4)
hold all;
plot(nv.T, nv.out('K_e'))
xlabel('Time [s]'); ylabel('K_e')
xlim([100 400])

subplot(3,2,5)
hold all;
plot(nv.T, nv.out('O2'))
xlabel('Time [s]'); ylabel('O2')
xlim([100 400])

subplot(3,2,6)
hold all;
plot(nv.T, nv.out('current'))
xlabel('Time [s]'); ylabel('current')
xlim([100 400])

figure(7);
subplot(2,2,1)
hold all;
plot(nv.T, nv.out('J_KIR_i'))
xlabel('Time [s]'); ylabel('KIR flux')
xlim([100 400])

subplot(2,2,2)
hold all;
plot(nv.T, nv.out('J_BK_k'))
xlabel('Time [s]'); ylabel('BK flux')
xlim([100 400])

subplot(2,2,3)
hold all;
plot(nv.T, nv.out('Ca_i'))
xlabel('Time [s]'); ylabel('Ca_i')
xlim([100 400])
