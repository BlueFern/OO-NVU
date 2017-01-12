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


odeopts = odeset('RelTol', 1e-04, 'AbsTol', 1e-04, 'MaxStep', 0.5, 'Vectorized', 1);

% Things to change:
XLIM1 = 1025; XLIM2 = 1200;      % x limits of any plots
PRESSURE = 4000;              % Pressure in Pa, change here!! Because the parameter is used in both SMCEC and WallMechanics ***********, 4000 Pa is default
PRESSURE_CHANGE=1500;         % Pressure at 100 sec till 300 seconds when switch is on, 2000 Pa is default
START_PULSE  =0;           % Start time of neuronal stimulation
LENGTH_PULSE = 2000;           % Length of stimulation 
GLU_SWITCH   = 1;             % Turn on glutamate input to SC, default 1
NO_SWITCH    = 1;             % Turn on Nitric Oxide production, default 1
TRPV4_SWITCH = 1;             % Turn on TRPV4 channel flux in AC, default 1
PLC_SWITCH   = 1;             % Turn on PLC change in time, default 1
PressureSwitch =0;            % Turn on pressure change in time, default 1
J_PLC_steadystate = 0.15;     % Jplc value in EC: 0.18 for steady state (default)
J_PLC_vasomotion  = 0.30;     % Jplc value in EC:  0.4 for vasomotion   (default) from 100 to 300 seconds

% Initiate the NVU
nv = NVU(Neuron('startpulse', START_PULSE, 'lengthpulse', LENGTH_PULSE, 'GluSwitch', GLU_SWITCH, 'NOswitch', NO_SWITCH), ...
    Astrocyte('startpulse', START_PULSE, 'lengthpulse', LENGTH_PULSE, 'rhoSwitch', GLU_SWITCH, 'TRPV4switch', TRPV4_SWITCH), ...
    WallMechanics('PressureSwitch', PressureSwitch,'P_T', PRESSURE,'Pressure_change', PRESSURE_CHANGE), ...
    SMCEC('PressureSwitch', PressureSwitch, 'J_PLC_steadystate', J_PLC_steadystate, 'J_PLC_vasomotion', J_PLC_vasomotion,  ...
    'NOswitch', NO_SWITCH, 'P_T', PRESSURE, 'PLCSwitch', PLC_SWITCH, 'Pressure_change', PRESSURE_CHANGE), 'odeopts', odeopts);

nv.T = linspace(0, XLIM2, 20000);     % Initiate time vector
nv.simulate();                       % Run simulation

%% Plot stuff
% Either run plot_all_variables.m to get everything, or write your own code
% here to plot specific variables or fluxes. 
% E.g.

figure(15);         % Initialise figure
data = ['R          ';'Ca_i       ';'I_i        ';'J_NaCa_i   ';'J_K_i      ';...
        'J_stretch_i';'J_KIR_i    '];
celldata = cellstr(data);
% factor=[-1 -1 -2 -1 -1 -1 -1]factor(i)*1970*
for i=[1:length(celldata)]
    subplot(3,3,i+1)
    hold all
    plot(nv.T, nv.out(char(celldata(i))))
    h(i)=gca();
    xlabel('Time [s]'); title(char(celldata(i)))
    xlim([XLIM1 XLIM2])
end

subplot(3,3,1)
hold all
plot(nv.T, nv.out('v_i'))
h(i)=gca();
xlabel('Time [s]'); title('v_i')
xlim([XLIM1 XLIM2])


subplot(3,3,9)
hold all
plot(nv.T, nv.out('V_coup_i'))
h(i)=gca();
xlabel('Time [s]'); title('V_coup_i')
xlim([XLIM1 XLIM2])

linkaxes(h, 'x');
