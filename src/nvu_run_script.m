%% Demonstration script for new NVU model

%% Version control
% This is the latest version of the NVU model: Version 2.1
% K+, NO, Astrocytic Ca2+, TRPV4, ECS, Neuron

%% Construct NVU
% The NVU consists of a number of submodules, implemented as MATLAB
% classes, presently a neuron/ECS, an astrocyte, a lumped SMC/EC model, and a model of
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

timeStart = datetime('now');
fprintf('Start time is %s\n', char(timeStart));

odeopts = odeset('RelTol', 1e-04, 'AbsTol', 1e-04, 'MaxStep', 0.5, 'Vectorized', 1);

%%%%%%%%%%%%%%%%%%%%%%% 
% Things you can change
%%%%%%%%%%%%%%%%%%%%%%%

FIG_NUM = 1;    % For plotting purposes
XLIM1 = 95;     
XLIM2 = 150;    % End of simulation

% Decide whether we want normal NVC (0) or CSD conditions (1)
CSD_SWITCH    = 0;  

CURRENT_TYPE  = 1;  % Types of current input:
                    % 1: normal block
                    % 2: two block stimulations (second stimulation is 8 sec after and 1 sec long)
                    % 3: obtained from whisker pad experimental input data (csv file)
                    % 4: whisker pad (from experiment) + locus coeruleus (pain pathway)

% For current type 1 or 2 use max current strength 0.022
% For current type 3 use max current strength 0.042
% For current type 4 use max current strength 0.035
% 
CURRENT_STRENGTH    = 0.022;    % Max strength of current input in mA/cm2

NEURONAL_START      = 100;      % Start of neuronal stimulation

% Used if CURRENT_STRENGTH = 1 or 2
NEURONAL_END        = 116;      % End of neuronal stimulation 

% Used if CURRENT_STRENGTH = 3 or 4
ISI = 7;                        % INDEX for time period between stimulations [0.6,1,2,3,4,6,8]
stim = 1;                       % INDEX for length of initial stimulation [2,8,16]

J_PLC           = 0.11;     % Jplc value in EC: 0.11 for steady state, 0.3 for oscillations
GLU_SWITCH      = 1;        % Turn on glutamate input (for NO and Ca2+ pathways)
NO_PROD_SWITCH  = 0;        % Turn on Nitric Oxide production 
TRPV_SWITCH     = 1;        % Turn on TRPV4 Ca2+ channel from AC to PVS
O2SWITCH        = 1;        % 0: ATP is plentiful, 1: ATP is limited (oxygen-limited regime, default)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Things you can change if you know what you're doing!
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Load initial NVU depending on whether we want normal NVC or CSD conditions

if CSD_SWITCH == 0
    nv = NVU(Neuron('k_syn', 11.5, 'CurrentType', CURRENT_TYPE, 'O2switch', O2SWITCH, 'startpulse', ...
        NEURONAL_START, 'lengthpulse', NEURONAL_END - NEURONAL_START, 'Istrength', CURRENT_STRENGTH, ...
        'GluSwitch', GLU_SWITCH, 'NOswitch', NO_PROD_SWITCH, ...
        'gNaleak_sa', 6.2378e-5, 'gKleak_sa', 2.1989e-4, 'gleak_sa', 10*6.2378e-5, 'gNaleak_d', 6.2961e-5, 'gKleak_d', 2.1987e-4, 'gleak_d', 10*6.2961e-5, 'Imax', 0.013*6), ...
        Astrocyte('trpv_switch', TRPV_SWITCH), ...
        WallMechanics(), ...
        SMCEC('J_PLC', J_PLC, 'NOswitch', NO_PROD_SWITCH), 'odeopts', odeopts);
else
    nv = NVU(Neuron('k_syn', 1, 'CurrentType', 1, 'O2switch', O2SWITCH, 'startpulse', ...
        NEURONAL_START, 'lengthpulse', 1, 'Istrength', 0.006, ...
        'GluSwitch', GLU_SWITCH, 'NOswitch', NO_PROD_SWITCH, ...
        'gNaleak_sa', 9.5999e-6, 'gKleak_sa', 3.4564e-5, 'gleak_sa', 10*9.5999e-6, 'gNaleak_d', 1.0187e-5, 'gKleak_d', 3.4564e-5, 'gleak_d', 10*1.0187e-5, 'Imax', 0.013), ...
        Astrocyte('trpv_switch', TRPV_SWITCH), ...
        WallMechanics(), ...
        SMCEC('J_PLC', J_PLC, 'NOswitch', NO_PROD_SWITCH), 'odeopts', odeopts);
end

% Adjust time vector
nv.neuron.params.dt = 0.001; dt = nv.neuron.params.dt;
nv.T = 0:dt:XLIM2;
numTimeSteps = length(nv.T);

%% Load whisker stimulation input data from file, save as I_Wh and put into Neuron.m
if nv.neuron.params.CurrentType == 3 || nv.neuron.params.CurrentType == 4
    
    load Zheng2010_data.mat
    actual_ISI = info.isi_duration(ISI);
    actual_stim = info.condition_stim_duration(stim);
    sum_neural_wh = zeros(size(neural_tim_vector));
    for animal = 1:11
        for experiment = 1:10
            sum_neural_wh = sum_neural_wh+neural_data(:,ISI,stim,experiment,animal)'; % Sum all data 
        end
    end
    mean_neural_wh = sum_neural_wh./110;  % Average the neural data over all animals and experiments, animals*experiments=110
    neural_tim_vector_shifted = neural_tim_vector + NEURONAL_START;    % Shift so stimulation begins at NEURONAL_START
    interp_neural_wh = interp1(neural_tim_vector_shifted, mean_neural_wh, nv.T); % Interpolate so there is data for all timesteps for NVU
    interp_neural_wh(isnan(interp_neural_wh))=0.02;   % Remove NaNs     
    
    nv.neuron.input_data = interp_neural_wh;   % Replace dummy in Neuron with input data, leave it at that for CURRENT_TYPE = 3
    I_Wh = interp_neural_wh;                   % Save as whisker pad current I_Wh
end

%% Construct additional pain pathway stimulation I_LC and input total current I_total to Neuron.m
if nv.neuron.params.CurrentType == 4
    
    % Construct first stimulation, set to work for any NEURONAL_START or any initial stimulus duration (2,8,16 s)
    time_1 = linspace(0, actual_stim, 10000); 
    I_1 = 0.0006 * time_1.^2 + 0.1;             % Chosen for the shape! Can be modified if needed
    time_1 = time_1 + NEURONAL_START;           % Shift so starts at NEURONAL_START
    I_1 = interp1(time_1, I_1, nv.T);           % Make same size as nv.T
    I_1(isnan(I_1)) = 0;                        % Replace NaNs with zeros
    
    % Construct second stimulation with duration 2 sec
    time_2 = linspace(0, 1, 10000);
    I_2 = 0.1*ones(size(time_2));
    time_2 = time_2 + NEURONAL_START + actual_stim + actual_ISI;
    I_2 = interp1(time_2, I_2, nv.T);
    I_2(isnan(I_2)) = 0;
    
    % Add together
    I_LC = I_1 + I_2;
    
    % Total current (whisker pad plus LC)
    I_total = I_Wh + I_LC;
    
    % Input to Neuron.m
    nv.neuron.input_data = I_total;
end

%% Run the simulation
nv.simulate() 

% % Run these 3 lines to take the ICs as the end of the last simulation run
% i.e. steady state ICs  
% ICs = (nv.U(end, :))';
% nv.u0 = ICs;
% nv.simulateManualICs() 

% Find a point in time 20 sec before neuronal stimulation has begun (preNeuronalStimTime1)
% and the point in time when stimulation begins (preNeuronalStimTime2) and get the values for 
% CBF, CBV, HBR, CMRO2 at that time, then find the midpoint between the min and max and normalise with that value. 
% Done in this way so that it also works when the variables are oscillatory (i.e. J_PLC = 0.3)
np = nv.neuron.params; % Shortcut for neuron parameters
preNeuronalStimTime1 = floor((NEURONAL_START-20)*numTimeSteps/XLIM2);
preNeuronalStimTime2 = floor((NEURONAL_START)*numTimeSteps/XLIM2);
CBF = nv.out('CBF'); CBV = (nv.out('CBV'))'; HBR = (nv.out('HbR'))'; CMRO2 = nv.out('CMRO2');
CBF_0 = 0.5*( max(CBF(preNeuronalStimTime1:preNeuronalStimTime2)) + min(CBF(preNeuronalStimTime1:preNeuronalStimTime2)) );
CBV_0 = 0.5*( max(CBV(preNeuronalStimTime1:preNeuronalStimTime2)) + min(CBV(preNeuronalStimTime1:preNeuronalStimTime2)) );
HBR_0 = 0.5*( max(HBR(preNeuronalStimTime1:preNeuronalStimTime2)) + min(HBR(preNeuronalStimTime1:preNeuronalStimTime2)) );
CMRO2_0 = 0.5*( max(CMRO2(preNeuronalStimTime1:preNeuronalStimTime2)) + min(CMRO2(preNeuronalStimTime1:preNeuronalStimTime2)) );

%Normalised variables
CBF_N = CBF./CBF_0;
CBV_N = CBV./CBV_0;
HBR_N = HBR./HBR_0;
CMRO2_N = CMRO2./CMRO2_0;

HBT_N = CBF_N .* HBR_N ./ CMRO2_N;                                          % Total hemoglobin (normalised)
HBO_N = (HBT_N - 1) - (HBR_N - 1) + 1;                                      % Oxyhemoglobin (normalised)
BOLD_N = 100 * np.V_0 * ( np.a_1 * (1 - HBR_N) - np.a_2 * (1 - CBV_N) );    % BOLD (percentage increase from 0)


%% Plot experimental and model CBF from data file
if nv.neuron.params.CurrentType == 3 || nv.neuron.params.CurrentType == 4 
    sum_cbf = zeros(size(cbf_tim_vector));
    for animal = 1:11
        for experiment = 1:10
            sum_cbf = sum_cbf+cbf_data(:,ISI,stim,experiment,animal)';
        end
    end
    mean_cbf = (sum_cbf./110) - 1;
    cbf_tim_vector_shifted = cbf_tim_vector + NEURONAL_START;    % Shift so stimulation begins at NEURONAL_START
    % Plot
    figure;
    plot(cbf_tim_vector_shifted, mean_cbf, ':k', nv.T, (nv.out('CBF')-CBF_0)./CBF_0, '-k', 'LineWidth', 1);
    ylabel('\Delta CBF')
    xlabel('Time [s]')
    xlim([90 150])
    ylim([-0.05 0.3])
    legend('experiment','model')
    title(['CBF with initial duration ' num2str(actual_stim) ', ISI ' num2str(actual_ISI)] );
    p1=patch([100 100+actual_stim 100+actual_stim 100],[-0.05 -0.05 0.3 0.3],'k');
    set(p1,'FaceAlpha',0.1,'EdgeColor', 'none');
    p2=patch([100+actual_stim+actual_ISI 100+actual_stim+actual_ISI+1 100+actual_stim+actual_ISI+1 100+actual_stim+actual_ISI],[-0.05 -0.05 0.3 0.3],'k');
    set(p2,'FaceAlpha',0.1,'EdgeColor', 'none');
end
%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot figures - choose your own variables, examples are listed below
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
figure(FIG_NUM);
subplot(4,2,1);
    hold all;
    plot(nv.T, nv.out('Ca_k'), 'LineWidth', 1);
    ylabel('Ca_k [\muM]');
    xlim([XLIM1 XLIM2])
subplot(4,2,2);
    hold all;
    plot(nv.T, nv.out('eet_k'), 'LineWidth', 1);
    ylabel('eet_k [\muM]');
    xlim([XLIM1 XLIM2])
subplot(4,2,3);
    hold all;
    plot(nv.T, nv.out('v_k'), 'LineWidth', 1);
    ylabel('v_k [mV]');
    xlim([XLIM1 XLIM2])
subplot(4,2,4);
    hold all;
    plot(nv.T, nv.out('w_k'), 'LineWidth', 1);
    ylabel('w_k [-]');
    xlim([XLIM1 XLIM2])
subplot(4,2,5);
    hold all;
    plot(nv.T, nv.out('K_p'), 'LineWidth', 1);
    ylabel('K_p [mM]');
    xlim([XLIM1 XLIM2])
subplot(4,2,6);
    hold all;
    plot(nv.T, nv.out('R'), 'LineWidth', 1);
    ylabel('Radius [\mum]');
    xlim([XLIM1 XLIM2])
subplot(4,2,7);
    hold all;
    plot(nv.T, nv.out('K_e'), 'LineWidth', 1);
    ylabel('K_e [mM]');
    xlim([XLIM1 XLIM2])
subplot(4,2,8);
    hold all;
    plot(nv.T, nv.out('K_s'), 'LineWidth', 1);
    ylabel('K_s [\muM]');
    xlim([XLIM1 XLIM2])
    

figure(FIG_NUM + 1);
subplot(1,2,1)
    plot(nv.T, nv.out('K_e'))
    xlim([90 130])
    ylim([3 7])
    xlabel('Time [s]')
    ylabel('K_e [mM]')
    p1=patch([100 116 116 100],[3 3 7 7],'k');
    set(p1,'FaceAlpha',0.1,'EdgeColor', 'none');  
subplot(1,2,2)
    plot(nv.T, nv.out('K_s')/1e3)
    xlim([90 130])
    ylim([0 15])
    xlabel('Time [s]')
    ylabel('K_s [mM]')
    p1=patch([100 116 116 100],[0 0 15 15],'k');
    set(p1,'FaceAlpha',0.1,'EdgeColor', 'none');
    
  
    
timeEnd = datetime('now');
fprintf('End time is %s\n', char(timeEnd));
