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

timeStart = datestr(now,'dd-mm-yyyy HH:MM:SS');
fprintf('Start time is %s\n', timeStart);

odeopts = odeset('RelTol', 1e-04, 'AbsTol', 1e-04, 'MaxStep', 0.5, 'Vectorized', 1);

XLIM1 = 290; XLIM2 = 380;
FIG_NUM = 3;


NEURONAL_START      = 300;      % Start of neuronal stimulation
NEURONAL_END        = 330;      % End of neuronal stimulation 
CURRENT_STRENGTH    = 0.015;  % Strength of current input in mA/cm2  0.015

ECS_START       = 100000000;      % Start of ECS K+ input
ECS_END         = 200000000;      % End of ECS K+ input
J_PLC           = 0.11;      % Jplc value in EC: 0.11 for steady state, 0.3 for oscillations

ECS             = 1;        % Include ECS compartment or not
CA_SWITCH       = 1;        % Turn on Ca2+ pathway
NO_INPUT_SWITCH = 1;        % Turn on NO stimulation
NO_PROD_SWITCH  = 1;        % Turn on Nitric Oxide production 
TRPV_SWITCH     = 1;        % Turn on TRPV4 Ca2+ channel from AC to PVS

nv = NVU(Neuron('v_switch', -50, 'startpulse', NEURONAL_START, 'lengthpulse', NEURONAL_END - NEURONAL_START, 'Istrength', CURRENT_STRENGTH, 'GluSwitch', NO_INPUT_SWITCH, 'NOswitch', NO_PROD_SWITCH, 't0_ECS', ECS_START), ...
    Astrocyte('trpv_switch', TRPV_SWITCH, 'startpulse', NEURONAL_START, 'lengthpulse', NEURONAL_END - NEURONAL_START, 't0_ECS', ECS_START, 'tend_ECS', ECS_END, 'GluSwitch', CA_SWITCH, 'ECSswitch', ECS, 'PVStoECS', 0, 'SCtoECS', 1), ...
    WallMechanics(), ...
    SMCEC('J_PLC', J_PLC, 'NOswitch', NO_PROD_SWITCH), 'odeopts', odeopts);

nv.T = linspace(0, XLIM2*2, 50000);
nv.simulate()

% % Run this to take the ICs as the end of the last simulation run i.e. steady state ICs
% ICs = (nv.U(end, :))';
% nv.u0 = ICs;
% nv.simulateManualICs() 

% Plot figures - whatever you want

% figure(2);
% plot(nv.T, nv.out('Glu'))
% ylabel('Glu')
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
    plot(nv.T, nv.out('CBV'), 'LineWidth', 1);
    ylabel('CBV');
    xlim([XLIM1 XLIM2])
subplot(3,3,7);
    hold all;
    plot(nv.T, nv.out('Glu'), 'LineWidth', 1);
    ylabel('Glutamate');
    xlim([XLIM1 XLIM2])
subplot(3,3,8);
    hold all;
    plot(nv.T, nv.out('E'), 'LineWidth', 1);
    ylabel('Oxygen extraction fraction');
    xlim([XLIM1 XLIM2])
subplot(3,3,9);
    hold all;
    plot(nv.T, nv.out('O2'), 'LineWidth', 1);
    ylabel('O2');
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

timeEnd = datestr(now,'dd-mm-yyyy HH:MM:SS');
fprintf('Finish time is %s\n', timeEnd);
