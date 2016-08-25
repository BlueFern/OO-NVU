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

clear all
clc
odeopts = odeset('RelTol', 1e-04, 'AbsTol', 1e-04, 'MaxStep', 0.5, 'Vectorized', 1);

XLIM1 = 998; XLIM2 = 1060;

FIG_NUM = 18;

ECS          = 1;        % Include ECS compartment or not
START_PULSE  = 1000;      % Time of first neuronal stimulation
LENGTH_PULSE = 1.5;       % Length of each stimulation 
J_PLC        = 0.18;     % Jplc value in EC: 0.18 for steady state, 0.4 for oscillations

nv = NVU(...
    Astrocyte('lengtht1', LENGTH_PULSE, 'F_input', 11.67, 'startpulse', START_PULSE, 'lengthpulse', LENGTH_PULSE, 'ECSswitch', ECS, 'PVStoECS', 0, 'SCtoECS', 1), ...
    WallMechanics(), ...
    SMCEC('J_PLC', J_PLC), 'odeopts', odeopts);

% Changes to input stimulus
nv.astrocyte.params.num_inputs = 100;     % Number of input stimuli
nv.astrocyte.params.input_spacing = 4;  % Time between each stimulus

% Changes to model to speed it up
nv.astrocyte.params.R_decay = 50*0.05;  % decay of K+ in PVS
nv.wall.params.wall_scaling = 10;  % speed up wall mechanics (old data based on pig, we model mice so probably way faster)



% nv.smcec.params.VR_sr = 14;  % based on Jai's model - cytosol 0.7pL : SR 0.05pS

% nv.smcec.params.C_i = 55;      % CICR rate
% nv.smcec.params.L_i = 0.025;      % SR leak rate
%nv.smcec.params.B_i = 2.025;    % SERCA pump rate
%nv.smcec.params.F_i = 0.23;     % IP3 channel strength



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
set(hAx(1), 'ylim', [19.3 19.8]);
%set(hAx(1), 'ylim', [16.5 20]);
%set(hAx(1), 'ylim', [20.3 20.8]);
set(hAx(2), 'ylim', [0 15]);
