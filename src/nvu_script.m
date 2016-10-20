%% Demonstration script for new NVU model version 1.2

%% Version control
% This is the latest version of the NVU model. It includes:
% TRPV4
% NO pathway
% Astrocytic Ca2+

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

% Things to change:
XLIM1 = 0; XLIM2 = 400;       % x limits of any plots
PRESSURE = 4000;              % Pressure in Pa, change here!! Because the parameter is used in both SMCEC and WallMechanics ***********, 4000 Pa is default
PRESSURE_CHANGE=2000;         % Pressure at 100 sec till 300 seconds when switch is on, 2000 Pa is default
START_PULSE  = 500;           % Start time of neuronal stimulation
LENGTH_PULSE = 100;           % Length of stimulation 
GLU_SWITCH   = 1;             % Turn on glutamate input to SC, default 1
NO_SWITCH    = 1;             % Turn on Nitric Oxide production, default 1
TRPV4_SWITCH = 1;             % Turn on TRPV4 channel flux in AC, default 1
PLC_SWITCH   = 0;             % Turn on PLC change in time, default 1
PressureSwitch = 1;           % Turn on pressure change in time, default 1
J_PLC_steadystate = 0.23;     % Jplc value in EC: 0.18 for steady state (default)
J_PLC_vasomotion  =  0.25;     % Jplc value in EC:  0.4 for vasomotion   (default) from 100 to 300 seconds


% Initiate the NVU
nv = NVU(Neuron('startpulse', START_PULSE, 'lengthpulse', LENGTH_PULSE, 'GluSwitch', GLU_SWITCH, 'NOswitch', NO_SWITCH), ...
    Astrocyte('startpulse', START_PULSE, 'lengthpulse', LENGTH_PULSE, 'rhoSwitch', GLU_SWITCH, 'TRPV4switch', TRPV4_SWITCH), ...
    WallMechanics('PressureSwitch', PressureSwitch,'P_T', PRESSURE,'Pressure_change', PRESSURE_CHANGE), ...
    SMCEC('PressureSwitch', PressureSwitch, 'J_PLC_steadystate', J_PLC_steadystate, 'J_PLC_vasomotion', ...
    J_PLC_vasomotion, 'NOswitch', NO_SWITCH, 'P_T', PRESSURE, 'PLCSwitch', PLC_SWITCH, 'Pressure_change', PRESSURE_CHANGE), 'odeopts', odeopts);

nv.T = linspace(0, XLIM2, 5000);     % Initiate time vector
nv.simulate();                       % Run simulation

%% Plot stuff
% Either run plot_all_variables.m to get everything, or write your own code
% here to plot specific variables or fluxes. 
% E.g.

figure(1);         % Initialise figure

subplot(4,4,1)
hold all
plot(nv.T, nv.out('I_j'))
h(1) = gca();
xlabel('Time [s]'); title('concentration IP3 EC')
xlim([XLIM1 XLIM2])
% line([108.9 108.9], ylim,'Color',[0,0.7,0.9]);
% line([121.9 121.9], ylim,'Color','r');
% line([119.5 119.5], ylim, 'Color', 'b');
% line([114.5 114.5], ylim, 'Color', 'g');

subplot(4,4,2)
hold all
plot(nv.T, nv.out('I_i'))
h(2) = gca();
xlabel('Time [s]'); title('concentration IP3 SMC')
xlim([XLIM1 XLIM2])
% line([108.9 108.9], ylim,'Color',[0,0.7,0.9]);
% line([121.9 121.9], ylim,'Color','r');
% line([119.5 119.5], ylim, 'Color', 'b');
% line([114.5 114.5], ylim, 'Color', 'g');

subplot(4,4,3)
hold all
plot(nv.T, nv.out('J_IP3_coup_i'))
h(3) = gca();
xlabel('Time [s]'); title('flux through IP3 coupling')
xlim([XLIM1 XLIM2])
% line([108.9 108.9], ylim,'Color',[0,0.7,0.9]);
% line([121.9 121.9], ylim,'Color','r');
% line([119.5 119.5], ylim, 'Color', 'b');
% line([114.5 114.5], ylim, 'Color', 'g');

subplot(4,4,4)
hold all
plot(nv.T, nv.out('J_IP3_i'))
h(4) = gca();
xlabel('Time [s]'); title('Ca2+ flux through IP3 activated SR')
xlim([XLIM1 XLIM2])
% line([108.9 108.9], ylim,'Color',[0,0.7,0.9]);
% line([121.9 121.9], ylim,'Color','r');
% line([119.5 119.5], ylim, 'Color', 'b');
% line([114.5 114.5], ylim, 'Color', 'g');

subplot(4,4,5)
hold all
plot(nv.T, nv.out('J_CICR_i'))
h(5) = gca();
xlabel('Time [s]'); title('Ca2+ flux through CICR channel in SMC [\muM/s]')
xlim([XLIM1 XLIM2])
% line([108.9 108.9], ylim,'Color',[0,0.7,0.9]);
% line([121.9 121.9], ylim,'Color','r');
% line([119.5 119.5], ylim, 'Color', 'b');
% line([114.5 114.5], ylim, 'Color', 'g');

subplot(4,4,6)
hold all
plot(nv.T, nv.out('Ca_i'))
h(6) = gca();
xlabel('Time [s]'); title('Ca2+ concentration SMC')
xlim([XLIM1 XLIM2])
% line([108.9 108.9], ylim,'Color',[0,0.7,0.9]);
% line([121.9 121.9], ylim,'Color','r');
% line([119.5 119.5], ylim, 'Color', 'b');
% line([114.5 114.5], ylim, 'Color', 'g');

subplot(4,4,7)
hold all
plot(nv.T, nv.out('F_r'))
h(7) = gca();
xlabel('Time [s]'); title('Fraction of attached myosin cross-bridges[-]')
xlim([XLIM1 XLIM2])
% line([108.9 108.9], ylim,'Color',[0,0.7,0.9]);
% line([121.9 121.9], ylim,'Color','r');
% line([119.5 119.5], ylim, 'Color', 'b');
% line([114.5 114.5], ylim, 'Color', 'g');

subplot(4,4,8)
hold all
plot(nv.T, nv.out('w_i'))
h(8) = gca();
xlabel('Time [s]'); title('Open probability of CA activated K channels')
xlim([XLIM1 XLIM2])
% line([108.9 108.9], ylim,'Color',[0,0.7,0.9]);
% line([121.9 121.9], ylim,'Color','r');
% line([119.5 119.5], ylim, 'Color', 'b');
% line([114.5 114.5], ylim, 'Color', 'g');

subplot(4,4,9)
hold all
plot(nv.T, nv.out('J_K_i'))
h(9) = gca();
xlabel('Time [s]'); title('K flux through SMC')
xlim([XLIM1 XLIM2])
% line([108.9 108.9], ylim,'Color',[0,0.7,0.9]);
% line([121.9 121.9], ylim,'Color','r');
% line([119.5 119.5], ylim, 'Color', 'b');
% line([114.5 114.5], ylim, 'Color', 'g');

subplot(4,4,10)
hold all
plot(nv.T, nv.out('v_i'))
h(10) = gca();
xlabel('Time [s]'); title('membrame potential SMC')
xlim([XLIM1 XLIM2])
% line([108.9 108.9], ylim,'Color',[0,0.7,0.9]);
% line([121.9 121.9], ylim,'Color','r');
% line([119.5 119.5], ylim, 'Color', 'b');
% line([114.5 114.5], ylim, 'Color', 'g');

subplot(4,4,11)
hold all
plot(nv.T, nv.out('J_VOCC_i'))
h(11) = gca();
xlabel('Time [s]'); title('Ca2+ flux through VOCC channel in SMC [\muM/s]')
xlim([XLIM1 XLIM2])
% line([108.9 108.9], ylim,'Color',[0,0.7,0.9]);
% line([121.9 121.9], ylim,'Color','r');
% line([119.5 119.5], ylim, 'Color', 'b');
% line([114.5 114.5], ylim, 'Color', 'g');

subplot(4,4,12)
hold all
plot(nv.T, nv.out('J_stretch_i'));
h(12) = gca();
xlabel('Time [s]'); title('Ca2+ flux through SMC stretch channel [\muM/s]')
xlim([XLIM1 XLIM2])
% line([108.9 108.9], ylim,'Color',[0,0.7,0.9]);
% line([121.9 121.9], ylim,'Color','r');
% line([119.5 119.5], ylim, 'Color', 'b');
% line([114.5 114.5], ylim, 'Color', 'g');

subplot(4,4,13)
hold all
plot(nv.T, nv.out('J_stretch_j'))
h(13) = gca();
xlabel('Time [s]'); title('Ca2+ flux through EC stretch channel [\muM/s]')
xlim([XLIM1 XLIM2])
% line([108.9 108.9], ylim,'Color',[0,0.7,0.9]);
% line([121.9 121.9], ylim,'Color','r');
% line([119.5 119.5], ylim, 'Color', 'b');
% line([114.5 114.5], ylim, 'Color', 'g');


subplot(4,4,14)
hold all
plot(nv.T, nv.out('s_i'))
h(14) = gca();
xlabel('Time [s]'); title('Ca2+ storage in SR')
xlim([XLIM1 XLIM2])
% line([108.9 108.9], ylim,'Color',[0,0.7,0.9]);
% line([121.9 121.9], ylim,'Color','r');
% line([119.5 119.5], ylim, 'Color', 'b');
% line([114.5 114.5], ylim, 'Color', 'g');

subplot(4,4,15)
hold all
plot(nv.T, nv.out('s_j'))
h(15) = gca();
xlabel('Time [s]'); title('Ca2+ storage in ER')
xlim([XLIM1 XLIM2])
% line([108.9 108.9], ylim,'Color',[0,0.7,0.9]);
% line([121.9 121.9], ylim,'Color','r');
% line([119.5 119.5], ylim, 'Color', 'b');
% line([114.5 114.5], ylim, 'Color', 'g');

subplot(4,4,16)
hold all
plot(nv.T, nv.out('R')*10^6)
h(16) = gca();
xlabel('Time [s]'); title('Radius of vessel [\mum]')
xlim([XLIM1 XLIM2])
% line([108.9 108.9], ylim,'Color',[0,0.7,0.9]);
% line([121.9 121.9], ylim,'Color','r');
% line([119.5 119.5], ylim, 'Color', 'b');
% line([114.5 114.5], ylim, 'Color', 'g');

linkaxes(h, 'x');



