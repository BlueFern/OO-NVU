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

XLIM1 = 100; XLIM2 = 600;
FIG_NUM = 18;

ECS          = 1;        % Include ECS compartment or not
START_PULSE  = 200;      % Time of neuronal stimulation
LENGTH_PULSE = 20;       % Length of stimulation 
J_PLC        = 0.18;     % Jplc value in EC: 0.18 for steady state, 0.4 for oscillations

nv = NVU(...
    Astrocyte('F_input', 2.67, 'startpulse', START_PULSE, 'lengthpulse', LENGTH_PULSE, 'ECSswitch', ECS, 'PVStoECS', 0, 'SCtoECS', 1), ...
    WallMechanics(), ...
    SMCEC('J_PLC', J_PLC), 'odeopts', odeopts);

nv.T = linspace(0, XLIM2, 5000);    

nv.simulate()


figure(FIG_NUM);

subplot(1,2,1)
hold all;
plot(nv.T, nv.out('K_p')/1e3);
xlabel('Time [s]'); title('K+ in PVS [mM]')
xlim([XLIM1 XLIM2])

subplot(1,2,2)
hold all;
plot(nv.T, nv.out('R')*1e6)
xlabel('Time [s]'); title('Radius [\mum]')
xlim([XLIM1 XLIM2])
