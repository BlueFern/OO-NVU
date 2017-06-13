%% Demonstration script for new NVU model

%% Version control
% This is the latest version of the NVU model: Version 2.0
% K+, NO, Astrocytic Ca2+, TRPV4, ECS, Neuron

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

XLIM1 = 270; XLIM2 = 360;
FIG_NUM = 1;


NEURONAL_START      = 300;      % Start of neuronal stimulation
NEURONAL_END        = 315;      % End of neuronal stimulation 
CURRENT_STRENGTH    = 0.025;  % Strength of current input in mA/cm2  (0.009 for subthreshold, 0.011 for bursting, 0.015 for constant stimulation)

ECS_START       = 30000;      % Start of ECS K+ input
ECS_END         = 32000;      % End of ECS K+ input

J_PLC           = 0.11;      % Jplc value in EC: 0.11 for steady state, 0.3 for oscillations

GLU_SWITCH      = 1;        % Turn on glutamate input (for NO and Ca2+ pathways)
NO_PROD_SWITCH  = 1;        % Turn on Nitric Oxide production 
TRPV_SWITCH     = 1;        % Turn on TRPV4 Ca2+ channel from AC to PVS
RK_SWITCH       = 0;        % Make R_k variable (1) or constant (0)
O2SWITCH        = 1;        % 0: ATP is plentiful, 1: ATP is limited (oxygen-limited regime, default)

% Output start time and estimated finish time based on how many sec of
% stimulation there will be and speed of computer 
%(home: 0.38, work: 0.88)
estimatedMinutes = (NEURONAL_END - NEURONAL_START)*0.88;
timeStart = datetime('now');
estimatedTimeEnd = timeStart + minutes(estimatedMinutes); 
fprintf('Start time is %s\nEstimated end time is %s\n', char(timeStart), char(estimatedTimeEnd));

nv = NVU(Neuron('SC_coup', 11.5, 'O2switch', O2SWITCH, 'startpulse', NEURONAL_START, 'lengthpulse', NEURONAL_END - NEURONAL_START, 'Istrength', CURRENT_STRENGTH, 'GluSwitch', GLU_SWITCH, 'NOswitch', NO_PROD_SWITCH, 't0_ECS', ECS_START, 'ECS_input', 9), ...
    Astrocyte('R_decay', 0.15, 'Rk_switch', RK_SWITCH, 'trpv_switch', TRPV_SWITCH, 'startpulse', NEURONAL_START, 'lengthpulse', NEURONAL_END - NEURONAL_START, 't0_ECS', ECS_START, 'tend_ECS', ECS_END), ...
    WallMechanics('wallMech', 1.1), ...
    SMCEC('J_PLC', J_PLC, 'NOswitch', NO_PROD_SWITCH), 'odeopts', odeopts);

numTimeSteps = 150000;
nv.T = linspace(0, XLIM2, numTimeSteps);
nv.simulate()

% Find a point in time 50 sec before neuronal stimulation has begun (preNeuronalStimTime1)
% and get the values for CBF, CBV, DHG at that time, then find the midpoint between 
% the min and max and normalise with that value. 
% Done in this way so that it also works when the variables are oscillatory (i.e. J_PLC = 0.3)
preNeuronalStimTime1 = floor((NEURONAL_START-50)*numTimeSteps/XLIM2);
preNeuronalStimTime2 = floor((NEURONAL_START)*numTimeSteps/XLIM2);
CBF = nv.out('CBF'); CBV = nv.out('CBV'); DHG = nv.out('DHG');
CBF_0 = 0.5*( max(CBF(preNeuronalStimTime1:preNeuronalStimTime2)) + min(CBF(preNeuronalStimTime1:preNeuronalStimTime2)) );
CBV_0 = 0.5*( max(CBV(preNeuronalStimTime1:preNeuronalStimTime2)) + min(CBV(preNeuronalStimTime1:preNeuronalStimTime2)) );
DHG_0 = 0.5*( max(DHG(preNeuronalStimTime1:preNeuronalStimTime2)) + min(DHG(preNeuronalStimTime1:preNeuronalStimTime2)) );

np = nv.neuron.params; % Shortcut for neuron parameters

% % Run this to take the ICs as the end of the last simulation run i.e. steady state ICs
% ICs = (nv.U(end, :))';
% nv.u0 = ICs;
% nv.simulateManualICs() 

% Plot figures - whatever you want
% % 
% figure(10);
% subplot(2,2,1);
% hold all;
% plot(t, state(:,ind.N_sa_Em))
% ylabel('v_{sa}')
% xlim([20 100])
% subplot(2,2,2);
% hold all;
% plot(t, (state(:,ind.N_K_s)/2.8e-8)/1000)
% ylabel('K_s')
% xlim([20 100])
% subplot(2,2,3);
% hold all;
% plot(t, state(:,ind.N_e_K))
% ylabel('K_e')
% xlim([20 100])
% subplot(2,2,4);
% hold all;
% plot(t, state(:,ind.R)*1e6)
% ylabel('R')
% xlim([20 100])

figure(FIG_NUM);
subplot(3,3,1);
    hold all;
    plot(nv.T, nv.out('v_sa'), 'LineWidth', 1);
    ylabel('v_{sa} [mV]');
    xlim([XLIM1 XLIM2])
subplot(3,3,2);
    hold all;
    plot(nv.T, nv.out('E'), 'LineWidth', 1);
    ylabel('Ozygen extraction fraction');
    xlim([XLIM1 XLIM2])
subplot(3,3,3);
    hold all;
    plot(nv.T, nv.out('O2'), 'LineWidth', 1);
    ylabel('O_2 [mM]');
    xlim([XLIM1 XLIM2])
subplot(3,3,4);
    hold all;
    plot(nv.T, (nv.out('CBF')-CBF_0)./CBF_0, 'LineWidth', 1);
    ylabel('\Delta CBF / CBF_0');
    xlim([XLIM1 XLIM2])
subplot(3,3,5);
    hold all;
    plot(nv.T, nv.out('R')*1e6, 'LineWidth', 1);
    ylabel('Radius [\mum]');
    xlim([XLIM1 XLIM2])
subplot(3,3,6);
    hold all;
    plot(nv.T, nv.out('CMRO2')./nv.out('CMRO2_init'), 'LineWidth', 1);
    ylabel('CMRO_2');
    xlim([XLIM1 XLIM2])
subplot(3,3,7);
    hold all;
    plot(nv.T, nv.out('CBV')./CBV_0, 'LineWidth', 1);
    ylabel('CBV');
    xlim([XLIM1 XLIM2])
subplot(3,3,8);
    hold all;
    plot(nv.T, nv.out('DHG')./DHG_0, 'LineWidth', 1);
    ylabel('DHG');
    xlim([XLIM1 XLIM2])
subplot(3,3,9);
    hold all;
    % BOLD using normalised variables
    plot(nv.T, 100 * np.V_0 * ( np.a_1 * (1 - nv.out('DHG')/DHG_0) - np.a_2 * (1 - nv.out('CBV')/CBV_0) ), 'LineWidth', 1);  
    ylabel('\Delta BOLD (%)');
    xlim([XLIM1 XLIM2])
    
%     % 
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
