%% Demonstration script for new NVU model

%% Version control
% This is the latest version of the NVU model: Version 1.2!
% K+, NO, Astrocytic Ca2+, TRPV4, ECS

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



odeopts = odeset('RelTol', 1e-04, 'AbsTol', 1e-04, 'MaxStep', 0.5, 'Vectorized', 1);

XLIM1 = 290; XLIM2 = 350;
FIG_NUM = 3;


NEURONAL_START      = 300;      % Start of neuronal stimulation
NEURONAL_END        = 310;      % End of neuronal stimulation 
CURRENT_STRENGTH    = 0.015;  % Strength of current input in mA/cm2  (0.009 for subthreshold, 0.011 for bursting, 0.015 for constant stimulation)

ECS_START       = 30000;      % Start of ECS K+ input
ECS_END         = 32000;      % End of ECS K+ input

J_PLC           = 0.11;      % Jplc value in EC: 0.11 for steady state, 0.35 for oscillations

GLU_SWITCH      = 1;        % Turn on glutamate input (for NO and Ca2+ pathways)
NO_PROD_SWITCH  = 1;        % Turn on Nitric Oxide production 
TRPV_SWITCH     = 1;        % Turn on TRPV4 Ca2+ channel from AC to PVS
RK_SWITCH       = 0;        % Make R_k variable (1) or constant (0)

% Output start time and estimated finish time based on how many sec of
% stimulation there will be and speed of computer 
%(home: 0.35, work: 0.85)
estimatedMinutes = (NEURONAL_END - NEURONAL_START)*0.35;
timeStart = datetime('now');
estimatedTimeEnd = timeStart + minutes(estimatedMinutes); 
fprintf('Start time is %s\nEstimated end time is %s\n', char(timeStart), char(estimatedTimeEnd));

nv = NVU(Neuron('SC_coup', 3, 'O2switch', 1, 'startpulse', NEURONAL_START, 'lengthpulse', NEURONAL_END - NEURONAL_START, 'Istrength', CURRENT_STRENGTH, 'GluSwitch', GLU_SWITCH, 'NOswitch', NO_PROD_SWITCH, 't0_ECS', ECS_START, 'ECS_input', 9), ...
    Astrocyte('Rk_switch', RK_SWITCH, 'trpv_switch', TRPV_SWITCH, 'startpulse', NEURONAL_START, 'lengthpulse', NEURONAL_END - NEURONAL_START, 't0_ECS', ECS_START, 'tend_ECS', ECS_END), ...
    WallMechanics('wallMech', 3), ...
    SMCEC('J_PLC', J_PLC, 'NOswitch', NO_PROD_SWITCH), 'odeopts', odeopts);

numTimeSteps = 150000;
nv.T = linspace(0, XLIM2, numTimeSteps);
nv.simulate()

% Find a point in time 10 sec before neuronal stimulation has begun and get
% the values for DHG and CBV at that time for normalisation of these two
% variables
preNeuronalStimTime = floor((NEURONAL_START-10)*numTimeSteps/XLIM2);
CBV_0 = nv.U(preNeuronalStimTime, nv.neuron.index.CBV);
DHG_0 = nv.U(preNeuronalStimTime, nv.neuron.index.DHG);

np = nv.neuron.params; % Shortcut for neuron parameters

% % Run this to take the ICs as the end of the last simulation run i.e. steady state ICs
% ICs = (nv.U(end, :))';
% nv.u0 = ICs;
% nv.simulateManualICs() 

% Plot figures - whatever you want
% % 
figure(20);
subplot(1,1,1);
hold all;
plot(nv.T, nv.out('K_e'), nv.T, nv.out('K_s')/1e3)
ylabel('Ke and Ks')
xlim([XLIM1 XLIM2])
% subplot(1,3,2);
% hold all;
% plot(nv.T, nv.out('rho'))
% ylabel('rho')
% xlim([XLIM1 XLIM2])
% subplot(1,3,3);
% hold all;
% plot(nv.T, nv.out('K_s')/1e3, [XLIM1 XLIM2],[6 6])
% ylabel('K_s')
% xlim([XLIM1 XLIM2])

figure(FIG_NUM);

subplot(3,3,1);
    hold all;
    plot(nv.T, nv.out('current'), 'LineWidth', 1);
    ylabel('current [mA/cm^2]');
    xlim([XLIM1 XLIM2])
subplot(3,3,2);
    hold all;
    plot(nv.T, nv.out('v_sa'), 'LineWidth', 1);
    ylabel('v_{sa} [mV]');
    xlim([XLIM1 XLIM2])
subplot(3,3,3);
    hold all;
    plot(nv.T, nv.out('K_e'), 'LineWidth', 1);
    ylabel('K_e [mM]');
    xlim([XLIM1 XLIM2])
subplot(3,3,4);
    hold all;
    plot(nv.T, nv.out('K_s')/1e3, 'LineWidth', 1);
    ylabel('K_s [mM]');
    xlim([XLIM1 XLIM2])
subplot(3,3,5);
    hold all;
    plot(nv.T, nv.out('R')*1e6, 'LineWidth', 1);
    ylabel('Radius [\mum]');
    xlim([XLIM1 XLIM2])
subplot(3,3,6);
    hold all;
    plot(nv.T, nv.out('CBV')/CBV_0, 'LineWidth', 1);
    ylabel('CBV');
    xlim([XLIM1 XLIM2])
subplot(3,3,7);
    hold all;
    plot(nv.T, nv.out('K_p'), 'LineWidth', 1);
    ylabel('K_p');
    xlim([XLIM1 XLIM2])
subplot(3,3,8);
    hold all;
    % BOLD using normalised variables
    plot(nv.T, 100 * np.V_0 * ( np.a_1 * (1 - nv.out('DHG')/DHG_0) - np.a_2 * (1 - nv.out('CBV')/CBV_0) ), 'LineWidth', 1);  
    ylabel('\Delta BOLD (%)');
    xlim([XLIM1 XLIM2])
subplot(3,3,9);
    hold all;
    plot(nv.T, nv.out('DHG')/DHG_0, 'LineWidth', 1);
    ylabel('Deoxyhemoglobin');
    xlim([XLIM1 XLIM2])
    
    % 
% figure(FIG_NUM + 1);
% 
% subplot(2,3,1)
%     hold all;
%     plot(nv.T, nv.out('Ca_k'), 'LineWidth', 1);
%     ylabel('Ca_k [\muM]');
%     xlim([XLIM1 XLIM2])
%     xlabel('Time [s]'); 
% subplot(2,3,2)
%     hold all;
%     plot(nv.T, nv.out('eet_k'), 'LineWidth', 1);
%     ylabel('eet_k [\muM]');
%     xlim([XLIM1 XLIM2])
%     xlabel('Time [s]'); 
% subplot(2,3,3)
%     hold all;
%     plot(nv.T, nv.out('w_k'), 'LineWidth', 1);
%     ylabel('w_k [-]');
%     xlim([XLIM1 XLIM2])
%     xlabel('Time [s]'); 
% subplot(2,3,4)
%     hold all;
%     plot(nv.T, nv.out('R')*1e6, 'LineWidth', 1);
%     ylabel('Radius [\mum]');
%     xlim([XLIM1 XLIM2])
%     xlabel('Time [s]'); 
% subplot(2,3,5)
%     hold all;
%     plot(nv.T, nv.out('K_p')/1e3, 'LineWidth', 1);
%     xlabel('Time [s]'); 
%     ylabel('K_p [mM]');
%     xlim([XLIM1 XLIM2])
% subplot(2,3,6)
%     hold all;
%     plot(nv.T, nv.out('R')*1e6, 'LineWidth', 1)
%     xlabel('Time [s]'); 
%     ylabel('Radius [\mum]');
%     xlim([XLIM1 XLIM2])

timeEnd = datetime('now');
fprintf('End time is %s\n', char(timeEnd));
