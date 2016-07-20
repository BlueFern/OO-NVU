%% Demonstration script for new NVU model

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
clear;

odeopts = odeset('RelTol', 1e-03, 'AbsTol', 1e-03, 'MaxStep', 1, 'Vectorized', 1);

nv = NVU(Astrocyte('startpulse', 200, 'lengthpulse', 20, 'PVStoECS', 0, 'SCtoECS', 1), ...
    WallMechanics(), ...
    SMCEC('J_PLC', 0.18), ...
    'odeopts', odeopts);

nv.T = linspace(0, 800, 1000);    

nv.simulate()


figure(9);
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
plot(nv.T, nv.out('K_p')/1e3)
xlabel('Time [s]'); ylabel('K_p')
