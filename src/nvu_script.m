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

nv = NVU(Neuron_test('startpulse', 100, 'lengthpulse', 15, 'GluSwitch', 0, 'NaTswitch', 0, 'buffSwitch', 0), ...
    Astrocyte_test('startpulse', 100, 'lengthpulse', 15, 'rhoSwitch', 1, 'blockSwitch', 1), ...
    WallMechanics(), ...
    SMCEC('J_PLC', 0.18), 'odeopts', odeopts);

nv.T = linspace(0, 300, 10000);    
    
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

subplot(3,2,2)
hold all;
plot(nv.T, nv.out('K_s')/1e3)
xlabel('Time [s]'); ylabel('K_s')

subplot(3,2,3)
hold all;
plot(nv.T, nv.out('v_sa'));
xlabel('Time [s]'); ylabel('v_{sa}')

subplot(3,2,4)
hold all;
plot(nv.T, nv.out('K_e'))
xlabel('Time [s]'); ylabel('K_e')

subplot(3,2,5)
hold all;
plot(nv.T, nv.out('O2'))
xlabel('Time [s]'); ylabel('O2')

subplot(3,2,6)
hold all;
plot(nv.T, nv.out('Ca_i'))
xlabel('Time [s]'); ylabel('Ca_i')
% 
% figure(10);
% plot(nv.T, nv.out('J_KIR_i'))
% xlabel('Time [s]'); ylabel('J_KIR_i')

% % Fix later to plot glutamate and rho
% subplot(4,2,6)
% hold all;
% plot(nv.T, nv.out('Glu'))
% xlabel('Time [s]'); ylabel('glutamate')
% 
% subplot(4,2,7)
% hold all;
% plot(nv.T, nv.out('rho'))
% xlabel('Time [s]'); ylabel('rho')
