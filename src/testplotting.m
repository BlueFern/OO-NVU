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

XLIM1 = 150; XLIM2 = 300;
FIG_NUM = 14;

START_PULSE  = 200;      % Time of neuronal stimulation
LENGTH_PULSE = 20;       % Length of stimulation 
GLU_SWITCH   = 1;        % Turn on glutamate input to SC
NO_SWITCH    = 1;        % Turn on Nitric Oxide production 
J_PLC        = 0.18;     % Jplc value in EC: 0.18 for steady state, 0.4 for oscillations

nv = NVU(Neuron('startpulse', START_PULSE, 'lengthpulse', LENGTH_PULSE, 'GluSwitch', GLU_SWITCH, 'NOswitch', NO_SWITCH), ...
    AstrocyteNoECS('startpulse', START_PULSE, 'lengthpulse', LENGTH_PULSE, 'rhoSwitch', GLU_SWITCH, 'PVStoECS', 0, 'SCtoECS', 1), ...
    WallMechanics(), ...
    SMCEC('J_PLC', J_PLC, 'NOswitch', NO_SWITCH), 'odeopts', odeopts);

nv.T = linspace(0, XLIM2, 5000);    

% BK conductance new, TRPV on
nv.astrocyte.params.r_buff = 0.05;
nv.astrocyte.params.G_BK_k = 225; % pS (later converted to mho m^-2)
nv.astrocyte.params.reverseBK = -0.08135; % Nernst potential of BK channel
nv.astrocyte.params.Ca_4 = 0.35; % uM
nv.astrocyte.params.v_7 = -13.57e-3; % V
nv.astrocyte.params.trpv_switch = 1;
nv.simulate()

figure(FIG_NUM);
subplot(3,2,1)
hold all;
plot(nv.T, nv.out('J_TRPV_k'))
xlabel('Time [s]'); title('TRPV Ca2+ flux [\muM/s]')
xlim([XLIM1 XLIM2])

subplot(3,2,2)
hold all;
plot(nv.T, nv.out('Ca_k'))
xlabel('Time [s]'); title('Ca2+ in Astrocyte [\muM]')
xlim([XLIM1 XLIM2])

subplot(3,2,3)
hold all;
plot(nv.T, nv.out('v_k'));
xlabel('Time [s]'); title('AC membrane potential [mV]')
xlim([XLIM1 XLIM2])

subplot(3,2,4)
hold all;
plot(nv.T, nv.out('w_k'))
xlabel('Time [s]'); title('BK Open Probability')
xlim([XLIM1 XLIM2])

subplot(3,2,5)
hold all;
plot(nv.T, nv.out('K_p'));
xlabel('Time [s]'); title('K+ in PVS [\muM]')
xlim([XLIM1 XLIM2])

subplot(3,2,6)
hold all;
plot(nv.T, nv.out('R'))
xlabel('Time [s]'); title('Radius [\mum]')
xlim([XLIM1 XLIM2])

%% BK conductance old, TRPV off
nv.astrocyte.params.switchBK = 1;
nv.astrocyte.params.reverseBK = 0; % V
nv.astrocyte.params.G_BK_k = 4300; % pS (later converted to mho m^-2)
nv.astrocyte.params.Ca_4 = 0.15; % uM
nv.astrocyte.params.v_7 = -15e-3; % V
nv.astrocyte.params.trpv_switch = 0;

nv.simulate()

figure(FIG_NUM);
subplot(3,2,1)
hold all;
plot(nv.T, nv.out('J_TRPV_k'))
xlabel('Time [s]'); title('TRPV Ca2+ flux [\muM/s]')
xlim([XLIM1 XLIM2])

subplot(3,2,2)
hold all;
plot(nv.T, nv.out('Ca_k'))
xlabel('Time [s]'); title('Ca2+ in Astrocyte [\muM]')
xlim([XLIM1 XLIM2])

subplot(3,2,3)
hold all;
plot(nv.T, nv.out('v_k'));
xlabel('Time [s]'); title('AC membrane potential [mV]')
xlim([XLIM1 XLIM2])

subplot(3,2,4)
hold all;
plot(nv.T, nv.out('w_k'))
xlabel('Time [s]'); title('BK Open Probability')
xlim([XLIM1 XLIM2])

subplot(3,2,5)
hold all;
plot(nv.T, nv.out('K_p'));
xlabel('Time [s]'); title('K+ in PVS [\muM]')
xlim([XLIM1 XLIM2])

subplot(3,2,6)
hold all;
plot(nv.T, nv.out('R'))
xlabel('Time [s]'); title('Radius [\mum]')
xlim([XLIM1 XLIM2])

%% 
legend('TRPV on, BK new','TRPV off, BK old')
