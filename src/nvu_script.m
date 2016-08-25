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

XLIM1 = 0; XLIM2 = 300;

FIG_NUM = 18;

ECS          = 1;        % Include ECS compartment or not
START_PULSE  = 100;      % Time of neuronal stimulation
LENGTH_PULSE = 100;       % Length of stimulation 
J_PLC        = 0.18;     % Jplc value in EC: 0.18 for steady state, 0.4 for oscillations

nv = NVU(...
    Astrocyte('lengtht1', 10, 'F_input', 2.67, 'startpulse', START_PULSE, 'lengthpulse', LENGTH_PULSE, 'ECSswitch', ECS, 'PVStoECS', 0, 'SCtoECS', 1), ...
    WallMechanics(), ...
    SMCEC('J_PLC', J_PLC), 'odeopts', odeopts);

nv.T = linspace(0, XLIM2, 50000);    
nv.simulate()

%% 
figure(FIG_NUM);
subplot(2,2,1)
hold all;
plot(nv.T, nv.out('K_s')/1e3);
xlabel('Time [s]'); title('K+ in SC [mM]')
xlim([XLIM1 XLIM2])

subplot(2,2,2)
hold all;
%rectangle('Position',[START_PULSE,18,LENGTH_PULSE,6],'FaceColor',[0.9,0.9,0.9],'EdgeColor','none')
plot(nv.T, nv.out('R')*1e6)
xlabel('Time [s]'); title('Radius [\mum]')
xlim([XLIM1 XLIM2])

subplot(2,2,3)
hold all;
plot(nv.T, nv.out('K_p')/1e3);
xlabel('Time [s]'); title('PVS K+ [mM]')
xlim([XLIM1 XLIM2])

subplot(2,2,4)
hold all;
plot(nv.T, nv.out('Ca_i'));
xlabel('Time [s]'); title('SMC Ca2+ [\muM]')
% plot(nv.T, nv.out('ft'));
% xlabel('Time [s]'); title('input')
xlim([XLIM1 XLIM2])

%% 
figure(FIG_NUM+11);
[hAx,hLine1,hLine2] = plotyy(nv.T, nv.out('R')*1e6, nv.T, nv.out('K_s')/1e3);
xlabel('Time [s]')
ylabel(hAx(1), 'Radius [um]')
ylabel(hAx(2), 'K+ in SC [mM]')
set(hAx, 'xlim', [1020, 1060]);
%set(hAx(1), 'ylim', [19.3 19.8]);
%set(hAx(1), 'ylim', [16.5 20]);
set(hAx(1), 'ylim', [20.3 20.8]);
set(hAx(2), 'ylim', [0 15]);


%%
figure(FIG_NUM+100);
subplot(2,2,1)
hold all;
plot(nv.T, nv.out('R')*1e6)
xlabel('Time [s]'); title('Radius [\mum]')
xlim([XLIM1 XLIM2])

subplot(2,2,2)
hold all;
plot(nv.T, nv.out('Ca_i'));
xlabel('Time [s]'); title('SMC Ca2+ [\muM]')
xlim([XLIM1 XLIM2])

subplot(2,2,3)
hold all;
plot(nv.T, nv.out('s_i'));
xlabel('Time [s]'); title('SR Ca2+ [\muM]')
xlim([XLIM1 XLIM2])



