%% Demonstration script for new NVU model

%% Version control
% This is the latest version of the NVU model. It includes:
% - basic K+ / Ca2+ 

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

nv = NVU(Neuron('startpulse', 200, 'lengthpulse', 200, 'KSwitch', 1, 'GluSwitch', 1), ...
    Astrocyte('startpulse', 200, 'lengthpulse', 200, 'rhoSwitch', 1, 'blockSwitch', 1, 'PVStoECS', 0, 'SCtoECS', 1), ...
    WallMechanics(), ...
    SMCEC('J_PLC', 0.18), 'odeopts', odeopts);

nv.T = linspace(0, 600, 5000);    
    
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

figure(7);
subplot(2,2,1)
hold all;
plot(nv.T, nv.out('R'))
xlabel('Time [s]'); ylabel('R')

subplot(2,2,2)
hold all;
plot(nv.T, nv.out('K_s')/1e3)
xlabel('Time [s]'); ylabel('K_s')

subplot(2,2,3)
hold all;
plot(nv.T, nv.out('Ca_k'));
xlabel('Time [s]'); ylabel('Ca_k')

subplot(2,2,4)
hold all;
plot(nv.T, nv.out('K_p')/1e3)
xlabel('Time [s]'); ylabel('K_p')

% subplot(3,2,5)
% hold all;
% plot(nv.T, nv.out('K_e')/1e3)
% xlabel('Time [s]'); ylabel('K_e')

% subplot(3,2,5)
% hold all;
% plot(nv.T, nv.out('O2'))
% xlabel('Time [s]'); ylabel('O2')

% 
% figure(5);
% hold on
% plot(nv.T, nv.out('K_e'))
% xlabel('time (s)')
% ylabel('K_e')
% hold off

% figure(6);
% hold on
% plot(nv.T, nv.out('J_lumen'))
% xlabel('time (s)')
% ylabel('J_{lumen}')
% hold off
% 
% figure(7);
% hold on
% plot(nv.T, nv.out('NO_j'))
% xlabel('time (s)')
% ylabel('NO_j')
% hold off
% 
% figure(8);
% hold on
% plot(nv.T, nv.out('p_NO_j'))
% xlabel('time (s)')
% ylabel('p_{NO_j}')
% hold off
% 
% figure(9);
% hold on
% plot(nv.T, nv.out('c_NO_j'))
% xlabel('time (s)')
% ylabel('c_{NO_j}')
% hold off
% 
% figure(10);
% hold on
% plot(nv.T, nv.out('p_NO_j') + nv.out('d_NO_j'))
% xlabel('time (s)')
% ylabel('test')
% hold off


