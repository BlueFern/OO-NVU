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
odeopts = odeset('RelTol', 1e-06, 'AbsTol', 1e-06, 'MaxStep', 1, 'Vectorized', 1);

% Handy parameters to modify
START_PULSE  = 200;         % Start time of neuronal stimulation
LENGTH_PULSE = 20;          % Length of stimulation
CURRENT_STRENGTH = 0.012;   % Strength of current input
J_PLC        = 0.18;        % Jplc value in EC: 0.18 for steady state, 0.4 for oscillations
XLIM1 = 150; XLIM2 = 250;   % Plotting limits


nv = NVU(Neuron('startpulse', START_PULSE, 'lengthpulse', LENGTH_PULSE, 'Istrength', CURRENT_STRENGTH), ...
    Astrocyte(), ...
    WallMechanics(), ...
    SMCEC('J_PLC', J_PLC), 'odeopts', odeopts);

nv.T = [0 300];    
    
nv.simulate()

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

figure(70);
subplot(4,2,1)
hold all;
plot(nv.T, nv.out('R')*1e6)
xlabel('Time [s]'); ylabel('R [\mum]')
xlim([XLIM1 XLIM2]);

subplot(4,2,2)
hold all;
plot(nv.T, nv.out('K_s')/1e3)
xlabel('Time [s]'); ylabel('K_s')
xlim([XLIM1 XLIM2]);

subplot(4,2,3)
hold all;
plot(nv.T, nv.out('v_sa'));
xlabel('Time [s]'); ylabel('v_{sa}')
xlim([XLIM1 XLIM2]);

subplot(4,2,4)
hold all;
plot(nv.T, nv.out('K_e'))
xlabel('Time [s]'); ylabel('K_e')
xlim([XLIM1 XLIM2]);

subplot(4,2,5)
hold all;
plot(nv.T, nv.out('K_p')/1e3)
xlabel('Time [s]'); ylabel('K_p')
xlim([XLIM1 XLIM2]);

% subplot(4,2,6)
% hold all;
% plot(nv.T, nv.out('Cl_sa'))
% xlabel('Time [s]'); ylabel('Cl soma/axon')
% xlim([XLIM1 XLIM2]);

% subplot(4,2,7)
% hold all;
% plot(nv.T, nv.out('CBF'))
% xlabel('Time [s]'); ylabel('CBF')
% xlim([380 500])

subplot(4,2,7)
hold all;
plot(nv.T, nv.out('current'))
xlabel('Time [s]'); ylabel('current')
xlim([XLIM1 XLIM2]);

% % 
% figure(7);
% subplot(2,2,1)
% hold all;
% plot(nv.T, nv.out('J_KIR_i'))
% xlabel('Time [s]'); ylabel('KIR flux')
% 
% subplot(2,2,2)
% hold all;
% plot(nv.T, nv.out('J_BK_k'))
% xlabel('Time [s]'); ylabel('BK flux')
% 
% subplot(2,2,3)
% hold all;
% plot(nv.T, nv.out('Ca_i'))
% xlabel('Time [s]'); ylabel('Ca_i')
