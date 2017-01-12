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
% clear all
% clc

odeopts = odeset('RelTol', 1e-04, 'AbsTol', 1e-04, 'MaxStep', 0.5, 'Vectorized', 1);

XLIM1 = 0; XLIM2 = 300;
FIG_NUM = 20;

ECS          = 1;        % Include ECS compartment or not
START_PULSE  = 1000;      % Time of neuronal stimulation
LENGTH_PULSE = 1000;       % Length of stimulation 
GLU_SWITCH   = 1;        % Turn on glutamate input to SC
NO_SWITCH    = 1;        % Turn on Nitric Oxide production 
K_SWITCH     = 1;        % Turn on K+ input into SC
J_PLC        = 0.1;     % Jplc value in EC: 0.18 for steady state, 0.4 for oscillations
J_STRETCH_FIX = -1;       % 1: old (wrong) , -1: new, J_stretch is - instead

nv = NVU(Neuron('F_input', 2.67, 'startpulse', START_PULSE, 'lengthpulse', LENGTH_PULSE, 'GluSwitch', GLU_SWITCH, 'NOswitch', NO_SWITCH, 'KSwitch', K_SWITCH), ...
    Astrocyte('ECS_input', 5, 'startpulse', START_PULSE, 'lengthpulse', LENGTH_PULSE, 'rhoSwitch', GLU_SWITCH, 'ECSswitch', ECS, 'PVStoECS', 0, 'SCtoECS', 1), ...
    WallMechanics(), ...
    SMCEC('J_PLC', J_PLC, 'NOswitch', NO_SWITCH, 'stretchFix', J_STRETCH_FIX), 'odeopts', odeopts);

nv.T = linspace(0, XLIM2*2, 5000);
%nv.u0 = temp;

% BK conductance new
% nv.astrocyte.params.r_buff = 0.05;
% nv.astrocyte.params.G_BK_k = 225; % pS (later converted to mho m^-2)
% nv.astrocyte.params.reverseBK = -0.08135; % Nernst potential of BK channel
% nv.astrocyte.params.Ca_4 = 0.35; % uM
% nv.astrocyte.params.v_7 = -13.57e-3; % V

nv.simulate()

% Run this to take the ICs as the end of the last simulation run
ICs = (nv.U(end, :))';
nv.u0 = ICs;
nv.simulateManualICs()

% figure(FIG_NUM);
% 
% subplot(1,2,1)
% hold all;
% plot(nv.T, nv.out('R')*1e6)
% xlabel('Time [s]'); ylabel('Radius [\mum]')
% xlim([XLIM1 XLIM2])
% 
% subplot(1,2,2)
% hold all;
% plot(nv.T, nv.out('w_k'));
% xlabel('Time [s]'); ylabel('w_k [-]')
% xlim([XLIM1 XLIM2])
% 
% figure(FIG_NUM+1)
% hold all;
% plot(nv.T, nv.out('R')*1e6)
% xlabel('Time [s]'); ylabel('Radius [\mum]')
% xlim([XLIM1 XLIM2])
% 
% % 
figure(FIG_NUM+2)

subplot(1,3,3)
hold all;
plot(nv.T, nv.out('R')*1e6)
xlabel('Time [s]'); ylabel('Radius [\mum]')
xlim([XLIM1 XLIM2])
ylim([16 30])

subplot(1,3,2)
hold all;
plot(nv.T, nv.out('K_p')/1e3)
xlabel('Time [s]'); ylabel('K_p (mM)')
xlim([XLIM1 XLIM2])
ylim([0 35])

subplot(1,3,1)
hold all;
plot(nv.T, nv.out('K_e')/1e3)
xlabel('Time [s]'); ylabel('K_e (mM)')
xlim([XLIM1 XLIM2])
ylim([0 55])

% subplot(2,2,4)
% hold all;
% plot(nv.T, nv.out('Ca_n'))
% xlabel('Time [s]'); title('Ca_n')
% xlim([XLIM1 XLIM2])
% 
% subplot(2,3,6)
% hold all;
% plot(nv.T, nv.out('J_BK_k'))
% xlabel('Time [s]'); title('BK flux')
% xlim([XLIM1 XLIM2])
% 
% subplot(1,2,3)
% hold all;
% plot(nv.T, nv.out('K_s')/1e3);
% xlabel('Time [s]'); title('K+ in SC [mM]')
% xlim([XLIM1 XLIM2])
% 
