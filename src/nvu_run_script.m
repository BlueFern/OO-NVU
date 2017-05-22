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

XLIM1 = 290; XLIM2 = 450;
FIG_NUM = 3;


NEURONAL_START      = 300;      % Start of neuronal stimulation
NEURONAL_END        = 340;      % End of neuronal stimulation 
CURRENT_STRENGTH    = 0.015;  % Strength of current input in mA/cm2  (0.009 for subthreshold, 0.011 for bursting, 0.015 for constant stimulation)

ECS_START       = 100000000;      % Start of ECS K+ input
ECS_END         = 200000000;      % End of ECS K+ input

J_PLC           = 0.11;      % Jplc value in EC: 0.11 for steady state, 0.3 for oscillations

GLU_SWITCH      = 1;        % Turn on glutamate input (for NO and Ca2+ pathways)
NO_PROD_SWITCH  = 1;        % Turn on Nitric Oxide production 
TRPV_SWITCH     = 1;        % Turn on TRPV4 Ca2+ channel from AC to PVS

% Output start time and estimated finish time based on how many sec of
% stimulation there will be and speed of computer 
%(home: 0.3, work: 0.86)
estimatedMinutes = (NEURONAL_END - NEURONAL_START)*0.86;
timeStart = datetime('now');
estimatedTimeEnd = timeStart + minutes(estimatedMinutes); 
fprintf('Start time is %s\nEstimated end time is %s\n', char(timeStart), char(estimatedTimeEnd));



nv = NVU(Neuron('v_switch', -50, 'startpulse', NEURONAL_START, 'lengthpulse', NEURONAL_END - NEURONAL_START, 'Istrength', CURRENT_STRENGTH, 'GluSwitch', GLU_SWITCH, 'NOswitch', NO_PROD_SWITCH, 't0_ECS', ECS_START), ...
    Astrocyte('trpv_switch', TRPV_SWITCH, 'startpulse', NEURONAL_START, 'lengthpulse', NEURONAL_END - NEURONAL_START, 't0_ECS', ECS_START, 'tend_ECS', ECS_END), ...
    WallMechanics('wallMech', 1), ...
    SMCEC('J_PLC', J_PLC, 'NOswitch', NO_PROD_SWITCH), 'odeopts', odeopts);

nv.T = linspace(0, XLIM2*2, 150000);
nv.simulate()

% % Run this to take the ICs as the end of the last simulation run i.e. steady state ICs
% ICs = (nv.U(end, :))';
% nv.u0 = ICs;
% nv.simulateManualICs() 

% Plot figures - whatever you want

figure(20);
subplot(1,3,1);
hold all;
plot(nv.T, nv.out('Glu'))
ylabel('Glu')
xlim([XLIM1 XLIM2])
subplot(1,3,2);
hold all;
plot(nv.T, nv.out('rho'))
ylabel('rho')
xlim([XLIM1 XLIM2])
subplot(1,3,3);
hold all;
plot(nv.T, nv.out('K_s')/1e3, [XLIM1 XLIM2],[6 6])
ylabel('K_s')
xlim([XLIM1 XLIM2])

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
    plot(nv.T, nv.out('CBV'), 'LineWidth', 1);
    ylabel('CBV');
    xlim([XLIM1 XLIM2])
% subplot(3,3,7);
%     hold all;
%     plot(nv.T, nv.out('NO_n'), 'LineWidth', 1);
%     ylabel('NO_n');
%     xlim([XLIM1 XLIM2])
subplot(3,3,8);
    hold all;
    plot(nv.T, nv.out('BOLD'), 'LineWidth', 1);
    ylabel('BOLD (%)');
    xlim([XLIM1 XLIM2])
subplot(3,3,9);
    hold all;
    plot(nv.T, nv.out('DHG'), 'LineWidth', 1);
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
