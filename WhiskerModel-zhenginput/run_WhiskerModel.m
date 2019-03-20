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

odeopts = odeset('RelTol', 1e-06, 'AbsTol', 1e-06, 'MaxStep', 0.5, 'Vectorized', 1);
FIG_NUM = 1;                % For plotting purposes
XLIM1 = 95000;
XLIM2 = 150000;             % End of simulation

NEURONAL_START      = 100000;      % Start of neuronal stimulation
NEURONAL_END        = 116000;      % End of neuronal stimulation

WHISKER_PAD = 0;            % 0: Model with square input pulse 
                            % 1: Model with whisker barrel input approximation with frequency 5 Hz from Zheng 
                            % 2: Model with constant individual input pulses from NeuroScience book in Zheng paradigm I, takes about 20 min to run

                            
ELASTICITY = 0;             % 0: Standard values (E_pas = 66e3 [Pa], E_act = 233e3 [Pa])
                            % 1: Increased radius change values (E_pas = 1e3 [Pa], E_act = 359e3 [Pa])
                            
ISI = 7;                    % INDEX for time period between stimulations [0.6,1,2,3,4,6,8]
stim = 3;                   % INDEX for length of initial stimulation [2,8,16]

CSD_SWITCH      = 0;
J_PLC           = 0.11e-3;     % Jplc value in EC: 0.11 for steady state, 0.3 for oscillations
GLU_SWITCH      = 1;        % Turn on glutamate input (for NO and Ca2+ pathways)
NO_PROD_SWITCH  = 1;        % Turn on Nitric Oxide production 
TRPV_SWITCH     = 1;        % Turn on TRPV4 Ca2+ channel from AC to PVS
O2SWITCH        = 1;        % 0: ATP is plentiful, 1: ATP is limited (oxygen-limited regime, default)

if WHISKER_PAD == 0
    actual_ISI = 0;
    secondpulselength = 0;
elseif WHISKER_PAD == 1 || WHISKER_PAD == 2     % Zheng input data for paradigm I:
    actual_ISI = 8000;                             % 8 seconds between stimulation periods
    secondpulselength = 1000;                      % 1 second stimulation for second pulse
end
secondpulse = NEURONAL_END + actual_ISI;        % Start of second pulse

% Alpha, beta and rest values adjusted for the different inputs
if WHISKER_PAD == 0
    alpha_K_e = 3.2;
    beta_K_e = 4.2e-3;
    alpha_Na_sa = 4.23;
    beta_Na_sa = 0.39e-3;
    alpha_Na_d = -2.12;
    beta_Na_d = 0.75e-3;
    K_eRest = 3.5;
    Na_saRest = 9.37;
    Na_dRest = 9.42;
elseif WHISKER_PAD == 1 || WHISKER_PAD == 2
    alpha_K_e = 3.5;
    beta_K_e = 10e-3;
    alpha_Na_sa = 4.23;
    beta_Na_sa = 0.65e-3;
    alpha_Na_d = -2.12;
    beta_Na_d = 0.7e-3;
    K_eRest = 3.567;
    Na_saRest = 9.237;
    Na_dRest = 9.318;
end

% Elastic modulus of the cerebral artery can be combined with subpopulation
% inputs
if ELASTICITY == 0
    E_PAS = 66e3;
    E_ACT = 233e3;
    P_INPUT = 1.0;
    Q_INPUT = 0.0;
elseif ELASTICITY == 1
    E_PAS = 1e3;
    E_ACT = 359e3;
    P_INPUT = 1.0;
    Q_INPUT = 0.0;
end

if CSD_SWITCH == 0
    nv = NVU(Neuron('Whiskerpad', WHISKER_PAD, 'startpulse', NEURONAL_START, ...
        'lengthpulse', NEURONAL_END - NEURONAL_START, 'secondpulse', secondpulse, ...
        'secondlength', secondpulselength, 'Pinput', P_INPUT, 'Qinput', Q_INPUT, ...
        'O2switch', O2SWITCH, 'GluSwitch', GLU_SWITCH, 'NOswitch', NO_PROD_SWITCH, ...
        'alpha_K_e', alpha_K_e, 'beta_K_e', beta_K_e, 'alpha_Na_sa', alpha_Na_sa, ...
        'beta_Na_sa', beta_Na_sa, 'alpha_Na_d', alpha_Na_d, 'beta_Na_d', beta_Na_d, ...
        'K_eBase', K_eRest, 'Na_saBase', Na_saRest, 'Na_dBase', Na_dRest), ...
        Astrocyte('trpv_switch', TRPV_SWITCH), ...
        WallMechanics('E_passive', E_PAS, 'E_active', E_ACT), ...
        SMCEC('J_PLC', J_PLC, 'NOswitch', NO_PROD_SWITCH), 'odeopts', odeopts);
else
    nv = NVU(Neuron('k_syn', 1, 'CurrentType', 1, 'O2switch', O2SWITCH, 'startpulse', ...
        NEURONAL_START, 'lengthpulse', 1, 'Istrength', 0.006, ...
        'GluSwitch', GLU_SWITCH, 'NOswitch', NO_PROD_SWITCH, ...
        'gNaleak_sa', 9.5999e-6, 'gKleak_sa', 3.4564e-5, 'gleak_sa', 10*9.5999e-6, 'gNaleak_d', ...
        1.0187e-5, 'gKleak_d', 3.4564e-5, 'gleak_d', 10*1.0187e-5, 'Imax', 0.013), ...
        Astrocyte('trpv_switch', TRPV_SWITCH), ...
        WallMechanics(), ...
        SMCEC('J_PLC', J_PLC, 'NOswitch', NO_PROD_SWITCH), 'odeopts', odeopts);
end

% Adjust time vector
nv.neuron.params.dt = 1; dt = nv.neuron.params.dt;                      % omzetten naar milliseconden?
nv.T = 0:dt:XLIM2;                                                          % nv.T is 100.000 stappen, maar sec
numTimeSteps = length(nv.T);

%% Create input data and put into Neuron.m
if WHISKER_PAD == 1   
    % Parameters
    Hz = 5;                     % [Hz or s^-1]                              % [individualpulses s^-1]
    peak_time = 1;              % [ms]                                      % Thalamic input increase time
    pulse_time = 15;            % [ms]                                      % Thalamic individual pulse input duration
    amplitude = 0.4;            % [spikes/ms]                               % Thalamic input amplitude
    base_input = 0.05;          % [spikes/ms]                               % Base thalamic input
    
    % Expressions
    msperpulse = 1000 ./ Hz;                                                % [ms individualpulse^-1]
    index_peak_time = peak_time + 1;                                        % Individual pulse increase length
    down_time = pulse_time - index_peak_time + 1;                           % Individual pulse decrease length
    rest_time = msperpulse - (down_time + index_peak_time) + 2;             % Rest time individual pulse
    numberOfPeriods = (NEURONAL_END - NEURONAL_START) ./ msperpulse;        % Nr of individual pulses in neuronal stim
    start_time = NEURONAL_START;
    end_time = XLIM2 - NEURONAL_END + 1;
    
    % Make start and end time resting values
    startSignal = linspace(base_input, base_input, start_time);
    endSignal = linspace(base_input, base_input, end_time);
    
    % Construct one individual pulse cycle
    risingSignal = linspace(base_input, amplitude, index_peak_time);
    fallingSignal = linspace(amplitude, base_input, down_time);
    restSignal = linspace(base_input, base_input, rest_time);

    % Combine rising, falling and resting sections into an individual pulse
    oneCycle = [risingSignal, fallingSignal(2:end-1), restSignal];
    
    % Now replicate this cycle several (numberOfPeriods) times.
    waveform = repmat(oneCycle, [1 numberOfPeriods]);
    
    nrofper = 10;
    waveform2 = repmat(oneCycle, [1 nrofper]);
    
    resting = repmat(0.05,1, 8000);
    
    % And put the values behind the baseline values before neuronal start
    timeline = [startSignal, waveform, endSignal];
    timeline2 = [startSignal, waveform, resting, waveform2, endSignal];
    
    nv.neuron.input_data = timeline2;
    
elseif WHISKER_PAD == 2
    % Parameters
    peak_time = 1;              % [ms]                                      % Thalamic input increase time
    pulse_time = 15;            % [ms]                                      % Thalamic individual pulse input duration
    amplitude = 0.4;            % [spikes/ms]                               % Thalamic input amplitude
    base_input = 0.05;          % [spikes/ms]                               % Base thalamic input
    
    % Expressions
    index_peak_time = peak_time + 1;                                        % Individual pulse increase length
    down_time = pulse_time - index_peak_time + 1;                           % Individual pulse decrease length
    rest_time = 0;
    numberOfPeriods = 1143;                                                 % Nr of individual pulses in 16 second input pulse
    start_time = NEURONAL_START;
    end_time = XLIM2 - NEURONAL_END + 1;
    
    % Make start and end time resting values
    startSignal = linspace(base_input, base_input, start_time);
    endSignal = linspace(base_input, base_input, end_time);
    
    % Construct one individual pulse cycle
    risingSignal = linspace(base_input, amplitude, index_peak_time);
    fallingSignal = linspace(amplitude, base_input, down_time);
    restSignal = linspace(base_input, base_input, rest_time);

    % Combine rising, falling and resting sections into an individual pulse
    oneCycle = [risingSignal, fallingSignal(2:end-1), restSignal];
    
    % Now replicate this cycle several (numberOfPeriods) times.
    waveform = repmat(oneCycle, [1 numberOfPeriods]);   
    nrofper = 134;
    waveform2 = repmat(oneCycle, [1 nrofper]);
    resting = repmat(0.05,1, 8000);
    
    % And put the values behind the baseline values before neuronal start
    timeline = [startSignal, waveform, endSignal];
    timeline2 = [startSignal, waveform, resting, waveform2, endSignal];
    
    nv.neuron.input_data = timeline2;
end
    
%% Run the simulation
nv.simulate()

%% Plot variables for article (model fitting)
Mathias_Wh_nv = load('Dormanns_Whiskerpad1_withsecondpulse_nv.mat');
Mathias_nv = struct2array(Mathias_Wh_nv);

Ke16 = Mathias_nv.out('K_e');
Ks16 = Mathias_nv.out('K_s');
Nasa16 = Mathias_nv.out('Na_sa');
Nad16 = Mathias_nv.out('Na_d');
dKedt16 = Mathias_nv.out('J_K_NEtoSC')/11.5;
 
figure(4)
hold on
plot(nv.T, nv.out('E_t'))
hold on
plot(nv.T, nv.out('E_t')-nv.out('I_t'))
hold on
plot(nv.T, nv.out('I_t'))
hold off
legend('E','E - I','I')
ylabel('Fraction')
xlim([XLIM1 130000])
figure
plot(nv.T, nv.out('I_t'))
ylabel('Proportion of inhibitory cells firing')
xlim([XLIM1 130000])

figure(2)
subplot(2,2,1);
    hold on
    plot(nv.T, nv.out('K_e'), 'LineWidth', 1);
    hold on
    plot(nv.T, Ke16, 'LineWidth', 1);
    hold off     
    ylabel('[K_e^+] [mM]');
    set(gca,'XTickLabel',0:10:30)
    xlabel('Time [s]')
    xlim([XLIM1 130000])
subplot(2,2,2);
    hold on
    plot(nv.T, nv.out('K_s'));
    hold on
    plot(nv.T, Ks16)
    hold off    
    ylabel('[K_s^+] [mM]');
    set(gca,'YTickLabel',0:5:25)
    set(gca,'XTickLabel',0:10:30)
    xlabel('Time [s]')
    xlim([XLIM1 130000])
subplot(2,2,3);
    hold on
    plot(nv.T, nv.out('Na_sa'), 'LineWidth', 1);
    hold on
    plot(nv.T, Nasa16);
    hold off   
    ylabel('[Na_{sa}^+] [mM]');
    set(gca,'XTickLabel',0:10:30)
    xlabel('Time [s]')
    xlim([XLIM1 130000])
subplot(2,2,4);
    hold on
    plot(nv.T, nv.out('Na_d'), 'LineWidth', 1);
    hold on
    plot(nv.T, Nad16)
    hold off
    ylabel('[Na_d^+] [mM]');
    set(gca,'XTickLabel',0:10:30)
    xlabel('Time [s]')
    xlim([XLIM1 130000])
    
figure(3)
    hold on
    plot(nv.T, nv.out('dKedt')*1000, 'LineWidth', 1);
%     hold on
%     plot(nv.T, dKedt16)
    hold off
    ylabel('^{d[K_e^+]}/_{dt} [mM/s]');
    set(gca,'XTickLabel',-5:5:30)
    xlabel('Time [s]')
    xlim([XLIM1 130000])

%% BOLD

% % Run these 3 lines to take the ICs as the end of the last simulation run
% i.e. steady state ICs  
% ICs = (nv.U(end, :))';
% nv.u0 = ICs;
% nv.simulateManualICs()

% Find a point in time 20 sec before neuronal stimulation has begun (preNeuronalStimTime1)
% and the point in time when stimulation begins (preNeuronalStimTime2) and get the values for 
% CBF, CBV, HBR, CMRO2 at that time, then find the midpoint between the min and max and normalise 
% with that value. 
% Done in this way so that it also works when the variables are oscillatory (i.e. J_PLC = 0.3)
np = nv.neuron.params; % Shortcut for neuron parameters
preNeuronalStimTime1 = floor((NEURONAL_START-20)*numTimeSteps/XLIM2);
preNeuronalStimTime2 = floor((NEURONAL_START)*numTimeSteps/XLIM2);
CBF = nv.out('CBF'); CBV = (nv.out('CBV'))'; HBR = (nv.out('HbR'))'; CMRO2 = nv.out('CMRO2');
CBF_0 = 0.5*( max(CBF(preNeuronalStimTime1:preNeuronalStimTime2)) + ...
    min(CBF(preNeuronalStimTime1:preNeuronalStimTime2)) );
CBV_0 = 0.5*( max(CBV(preNeuronalStimTime1:preNeuronalStimTime2)) + ...
    min(CBV(preNeuronalStimTime1:preNeuronalStimTime2)) );
HBR_0 = 0.5*( max(HBR(preNeuronalStimTime1:preNeuronalStimTime2)) + ...
    min(HBR(preNeuronalStimTime1:preNeuronalStimTime2)) );
CMRO2_0 = 0.5*( max(CMRO2(preNeuronalStimTime1:preNeuronalStimTime2)) + ...
    min(CMRO2(preNeuronalStimTime1:preNeuronalStimTime2)) );

%Normalised variables
CBF_N = CBF./CBF_0;
CBV_N = CBV./CBV_0;
HBR_N = HBR./HBR_0;
CMRO2_N = CMRO2./CMRO2_0;

HBT_N = CBF_N .* HBR_N ./ CMRO2_N;                                      % Total hemoglobin(normalised)
HBO_N = (HBT_N - 1) - (HBR_N - 1) + 1;                                  % Oxyhemoglobin(normalised)
BOLD_N = 100 * np.V_0 * ( np.a_1 * (1 - HBR_N) - np.a_2 * (1 - CBV_N) );% BOLD(percentage increase from 0)

%% Plot experimental, Whisker Barrel model and Mathias et al. model CBF from Zheng data file
if WHISKER_PAD == 1 || WHISKER_PAD == 2
    Mathias_Wh_nv = load('Dormanns_Whiskerpad1_withsecondpulse_nv.mat');
    Mathias_Wh_nv = struct2array(Mathias_Wh_nv);
elseif WHISKER_PAD == 0
    Mathias_Wh_nv = load('Dormanns_Squarepulse_withsecondpulse_nv.mat');
    Mathias_Wh_nv = struct2array(Mathias_Wh_nv);
end

Mathias_Wh_np = Mathias_Wh_nv.neuron.params; % Shortcut for neuron parameters
preNeuronalStimTime1 = floor((NEURONAL_START-20)*numTimeSteps/XLIM2);
preNeuronalStimTime2 = floor((NEURONAL_START)*numTimeSteps/XLIM2);
MWh_CBF = Mathias_Wh_nv.out('CBF'); MWh_CBV = (Mathias_Wh_nv.out('CBV'))'; ...
    MWh_HBR = (Mathias_Wh_nv.out('HbR'))'; MWh_CMRO2 = Mathias_Wh_nv.out('CMRO2');
MWh_CBF_0 = 0.5*( max(MWh_CBF(preNeuronalStimTime1:preNeuronalStimTime2)) + ...
    min(MWh_CBF(preNeuronalStimTime1:preNeuronalStimTime2)) );
MWh_CBV_0 = 0.5*( max(MWh_CBV(preNeuronalStimTime1:preNeuronalStimTime2)) + ...
    min(MWh_CBV(preNeuronalStimTime1:preNeuronalStimTime2)) );
MWh_HBR_0 = 0.5*( max(MWh_HBR(preNeuronalStimTime1:preNeuronalStimTime2)) + ...
    min(MWh_HBR(preNeuronalStimTime1:preNeuronalStimTime2)) );
MWh_CMRO2_0 = 0.5*( max(MWh_CMRO2(preNeuronalStimTime1:preNeuronalStimTime2)) + ...
    min(MWh_CMRO2(preNeuronalStimTime1:preNeuronalStimTime2)) );

%Normalised variables
MWh_CBF_N = MWh_CBF./MWh_CBF_0;
MWh_CBV_N = MWh_CBV./MWh_CBV_0;
MWh_HBR_N = MWh_HBR./MWh_HBR_0;
MWh_CMRO2_N = MWh_CMRO2./MWh_CMRO2_0;

MWh_HBT_N = MWh_CBF_N .* MWh_HBR_N ./ MWh_CMRO2_N;                                  % Total hemoglobin (normalised)
MWh_HBO_N = (MWh_HBT_N - 1) - (MWh_HBR_N - 1) + 1;                                % Oxyhemoglobin (normalised)
MWh_BOLD_N = 100 * Mathias_Wh_np.V_0 * ( Mathias_Wh_np.a_1 * (1 - MWh_HBR_N) - ...
    Mathias_Wh_np.a_2 * (1 - MWh_CBV_N) );                                       % BOLD (percentage increase from 0)

load Zheng2010_data.mat
sum_cbf = zeros(size(cbf_tim_vector));
for animal = 1:11
    for experiment = 1:10
        sum_cbf = sum_cbf+cbf_data(:,ISI,stim,experiment,animal)';
    end
end
mean_cbf = (sum_cbf./110) - 1;
cbf_tim_vector_shifted = cbf_tim_vector * 1000 + NEURONAL_START; % Shift so stimulation 
% begins at NEURONAL_START

if WHISKER_PAD == 0
    % Plot
%     figure(1)
%     hold all
%     plot(cbf_tim_vector_shifted, mean_cbf, ':k', nv.T, ...
%         ((nv.out('CBF')-CBF_0)./CBF_0), '--k', Mathias_Wh_nv.T * 1000, ...
%         ((Mathias_Wh_nv.out('CBF')-MWh_CBF_0)./MWh_CBF_0), '-k','LineWidth', 1);
%     ylabel('\Delta CBF')
%     xlabel('Time [s]')
%     xlim([90000 150000])
%     set(gca,'XTickLabel',-10:10:50)
%     ylim([-0.05 0.45])
%     legend('Zheng experiment','Wilson & Cowan','Mathias et al.')
%     title(['CBF with initial duration ' num2str(NEURONAL_END - NEURONAL_START)...
%         's, ISI ' num2str(actual_ISI) 's'] );
%     p1=patch([100000 116000 116000 100000],[-0.05 -0.05 0.5 0.5],'k');
%     set(p1,'FaceAlpha',0.1,'EdgeColor', 'none');
    
elseif WHISKER_PAD == 1 || WHISKER_PAD == 2
    % Plot
    figure(100)
    hold all
    % WITH MATHIAS ET AL. RESULTS
    plot(cbf_tim_vector_shifted, mean_cbf, ':k')
    plot(nv.T, ((nv.out('CBF')-CBF_0)./CBF_0), '--k')
    plot(Mathias_Wh_nv.T * 1000, ((Mathias_Wh_nv.out('CBF')-MWh_CBF_0)./MWh_CBF_0), '-k','LineWidth', 1);
    % WITHOUT MATHIAS ET AL. RESULTS
%     plot(cbf_tim_vector_shifted, mean_cbf, ':k', nv.T, ...
%         ((nv.out('CBF')-CBF_0)./CBF_0), '--k', 'LineWidth', 1);
    ylabel('\Delta CBF')
    xlabel('Time [s]')
    xlim([90000 150000])
    set(gca,'XTickLabel',-10:10:50)
    ylim([-0.05 0.35])
    title(['CBF with initial duration ' num2str((NEURONAL_END - NEURONAL_START)/1000)...
        's, ISI ' num2str(actual_ISI/1000) 's']);
    p1=patch([100000 116000 116000 100000],[-0.05 -0.05 0.5 0.5],'k');
    set(p1,'FaceAlpha',0.1,'EdgeColor', 'none');
    p2=patch([124000 125000 125000 124000],[-0.05 -0.05 0.5 0.5],'k');
    set(p2,'FaceAlpha',0.1,'EdgeColor', 'none');
    legend('Zheng experiment','Whisker model','Mathias et al.')
end

%% Plot experimental data by Berwick
if WHISKER_PAD == 1 || WHISKER_PAD == 2
    Mathias_nv = load('Dormanns_Whiskerpad1_withoutsecondpulse_nv.mat');
    Mathias_nv = struct2array(Mathias_nv);
elseif WHISKER_PAD == 0
    Mathias_nv = load('Dormanns_Squarpulse_withoutsecondpulse_nv.mat');
    Mathias_nv = struct2array(Mathias_nv);
end

Mathias_np = Mathias_nv.neuron.params; % Shortcut for neuron parameters
preNeuronalStimTime1 = floor((NEURONAL_START-20)*numTimeSteps/XLIM2);
preNeuronalStimTime2 = floor((NEURONAL_START)*numTimeSteps/XLIM2);
M_CBF = Mathias_nv.out('CBF'); M_CBV = (Mathias_nv.out('CBV'))'; ...
    M_HBR = (Mathias_nv.out('HbR'))'; M_CMRO2 = Mathias_nv.out('CMRO2');
M_CBF_0 = 0.5*( max(M_CBF(preNeuronalStimTime1:preNeuronalStimTime2)) + ...
    min(M_CBF(preNeuronalStimTime1:preNeuronalStimTime2)) );
M_CBV_0 = 0.5*( max(M_CBV(preNeuronalStimTime1:preNeuronalStimTime2)) + ...
    min(M_CBV(preNeuronalStimTime1:preNeuronalStimTime2)) );
M_HBR_0 = 0.5*( max(M_HBR(preNeuronalStimTime1:preNeuronalStimTime2)) + ...
    min(M_HBR(preNeuronalStimTime1:preNeuronalStimTime2)) );
M_CMRO2_0 = 0.5*( max(M_CMRO2(preNeuronalStimTime1:preNeuronalStimTime2)) + ...
    min(M_CMRO2(preNeuronalStimTime1:preNeuronalStimTime2)) );

%Normalised variables
M_CBF_N = M_CBF./M_CBF_0;
M_CBV_N = M_CBV./M_CBV_0;
M_HBR_N = M_HBR./M_HBR_0;
M_CMRO2_N = M_CMRO2./M_CMRO2_0;

M_HBT_N = M_CBF_N .* M_HBR_N ./ M_CMRO2_N;                                  % Total hemoglobin (normalised)
M_HBO_N = (M_HBT_N - 1) - (M_HBR_N - 1) + 1;                                % Oxyhemoglobin (normalised)
M_BOLD_N = 100 * Mathias_np.V_0 * ( Mathias_np.a_1 * (1 - M_HBR_N) - ...
    Mathias_np.a_2 * (1 - M_CBV_N) );                                       % BOLD (percentage increase from 0)

load scaledBerwickData.mat
scaled_time_oxy = NEURONAL_START + scaled_time_oxy * 1000;
scaled_time_deoxy = NEURONAL_START + scaled_time_deoxy * 1000;
scaled_time_hbt = NEURONAL_START + scaled_time_hbt * 1000;
figure(5)
hold all
% WITH MATHIAS
% plot(scaled_time_oxy, scaled_oxy, ':r', nv.T(95000:130000,1), HBO_N(1,95000:130000), '--r', ...
%     Mathias_nv.T(95000:130000,1)*1000, M_HBO_N(1,95000:130000), 'r', 'LineWidth', 1)
% plot(scaled_time_deoxy, scaled_deoxy, ':b', nv.T(95000:130000,1), HBR_N(1,95000:130000), '--b', ...
%     Mathias_nv.T(95000:130000,1)*1000, M_HBR_N(1,95000:130000), 'b', 'LineWidth', 1)
% plot(scaled_time_hbt, scaled_totaloxy, ':g', nv.T(95000:130000,1), HBT_N(1,95000:130000), '--g',...
%     Mathias_nv.T(95000:130000,1)*1000, M_HBT_N(1,95000:130000), 'g', 'LineWidth', 1)
% WITHOUT MATHIAS
% plot(scaled_time_oxy, scaled_oxy, ':r', nv.T(95000:130000,1), HBO_N(1,95000:130000), '--r', 'LineWidth', 1)
% plot(scaled_time_deoxy, scaled_deoxy, ':b', nv.T(95000:130000,1), HBR_N(1,95000:130000), '--b', 'LineWidth', 1)
% plot(scaled_time_hbt, scaled_totaloxy, ':g', nv.T(95000:130000,1), HBT_N(1,95000:130000), '--g', 'LineWidth', 1)
% xlim([95000 130000])
% set(gca,'XTickLabel',-5:5:30)
% p1=patch([100000 116000 116000 100000],[0.7 0.7 1.3 1.3],'k');
% set(p1,'FaceAlpha',0.05,'EdgeColor', 'none');
% ylim([0.94,1.1])
%ylim([0.98,1.042])
%legend('HbO experiment','HbO model','HbO Mathias model','HbR experiment',...
%    'HbR model','HbR Mathias model','HbT experiment','HbT model','HbT Mathias model')


figure;
hold all;
plot(scaled_time_oxy, scaled_oxy, 'r', 'LineWidth', 1)
plot(scaled_time_deoxy, scaled_deoxy, 'b', 'LineWidth', 1)
plot(scaled_time_hbt, scaled_totaloxy, 'g', 'LineWidth', 1)
xlim([95000 130000])
xlabel('Time [s]')
ylim([0.98 1.045]);
set(gca,'XTickLabel',-5:5:30)
p1=patch([100000 116000 116000 100000],[0.7 0.7 1.3 1.3],'k');
set(p1,'FaceAlpha',0.05,'EdgeColor', 'none');
legend('HbO','HbR','HbT');
title('Hemodynamics: Berwick Experiment')

figure;
hold all;
plot(nv.T(95000:130000,1), HBO_N(1,95000:130000), 'r', 'LineWidth', 1)
plot(nv.T(95000:130000,1), HBR_N(1,95000:130000), 'b', 'LineWidth', 1)
plot(nv.T(95000:130000,1), HBT_N(1,95000:130000), 'g', 'LineWidth', 1)
xlim([95000 130000])
ylim([0.98 1.01]);
xlabel('Time [s]')
set(gca,'XTickLabel',-5:5:30)
p1=patch([100000 116000 116000 100000],[0.7 0.7 1.3 1.3],'k');
set(p1,'FaceAlpha',0.05,'EdgeColor', 'none');
legend('HbO','HbR','HbT');
title('Hemodynamics: Whisker Model')

figure;
hold all;
plot(Mathias_nv.T(95000:130000,1)*1000, M_HBO_N(1,95000:130000), 'r', 'LineWidth', 1)
plot(Mathias_nv.T(95000:130000,1)*1000, M_HBR_N(1,95000:130000), 'b', 'LineWidth', 1)
plot(Mathias_nv.T(95000:130000,1)*1000, M_HBT_N(1,95000:130000), 'g', 'LineWidth', 1)
xlim([95000 130000])
ylim([0.9 1.17]);
xlabel('Time [s]')
set(gca,'XTickLabel',-5:5:30)
p1=patch([100000 116000 116000 100000],[0.7 0.7 1.3 1.3],'k');
set(p1,'FaceAlpha',0.05,'EdgeColor', 'none');
legend('HbO','HbR','HbT');
title('Hemodynamics: Mathias et al. Model')


%% Plot neuron variables for report
% Choose between square input pulses or whiskerpad from old nvu model

% Mathias_Wh_nv = load('Dormanns_Whiskerpad1_withsecondpulse_nv.mat');
% Mathias_Wh_nv = struct2array(Mathias_Wh_nv);
% Ke16 = Mathias_Wh_nv.out('K_e');
% Ks16 = Mathias_Wh_nv.out('K_s');
% Nasa16 = Mathias_Wh_nv.out('Na_sa');
% Nad16 = Mathias_Wh_nv.out('Na_d');
% dKedt16 = Mathias_Wh_nv.out('J_K_NEtoSC')/11.5;

% Mathias_Sq_nv = load('Dormanns_Squarpulse_withoutsecondpulse_nv.mat');
% Mathias_Sq_nv = struct2array(Mathias_Sq_nv);
% Ke16 = Mathias_Sq_nv.out('K_e');
% Ks16 = Mathias_Sq_nv.out('K_s');
% Nasa16 = Mathias_Sq_nv.out('Na_sa');
% Nad16 = Mathias_Sq_nv.out('Na_d');
% dKedt16 = Mathias_Sq_nv.out('J_K_NEtoSC')/11.5;
 
% figure(4)
% hold on
% plot(nv.T, nv.out('E_t'))
% hold on
% plot(nv.T, nv.out('E_t')-nv.out('I_t'))
% hold on
% plot(nv.T, nv.out('I_t'))
% hold off
% legend('E','E - I','I')
% ylabel('Fraction')
% xlim([XLIM1 XLIM2])
% figure
% plot(nv.T, nv.out('I_t'))
% ylabel('Proportion of inhibitory cells firing')
% xlim([XLIM1 XLIM2])
% 
% figure(2)
% subplot(2,2,1);
%     hold on
%     plot(nv.T, nv.out('K_e'), 'LineWidth', 1);
%     hold on
%     plot(nv.T, Ke16, 'LineWidth', 1);
%     hold off     
%     ylabel('[K_e^+] [mM]');
%     set(gca,'XTickLabel',0:10:30)
%     xlabel('Time [s]')
%     xlim([XLIM1 XLIM2])
% subplot(2,2,2);
%     hold on
%     plot(nv.T, nv.out('K_s'));
%     hold on
%     plot(nv.T, Ks16)
%     hold off    
%     ylabel('[K_s^+] [mM]');
%     set(gca,'YTickLabel',0:5:25)
%     set(gca,'XTickLabel',0:10:30)
%     xlabel('Time [s]')
%     xlim([XLIM1 XLIM2])
% subplot(2,2,3);
%     hold on
%     plot(nv.T, nv.out('Na_sa'), 'LineWidth', 1);
%     hold on
%     plot(nv.T, Nasa16);
%     hold off   
%     ylabel('[Na_{sa}^+] [mM]');
%     set(gca,'XTickLabel',0:10:30)
%     xlabel('Time [s]')
%     xlim([XLIM1 XLIM2])
% subplot(2,2,4);
%     hold on
%     plot(nv.T, nv.out('Na_d'), 'LineWidth', 1);
%     hold on
%     plot(nv.T, Nad16)
%     hold off
%     ylabel('[Na_d^+] [mM]');
%     set(gca,'XTickLabel',0:10:30)
%     xlabel('Time [s]')
%     xlim([XLIM1 XLIM2])
%     
% figure(3)
%     hold on
%     plot(nv.T, nv.out('dKedt')*1000, 'LineWidth', 1);
% %     hold on
% %     plot(nv.T, dKedt16)
%     hold off
%     ylabel('^{d[K_e^+]}/_{dt} [mM/s]');
%     set(gca,'XTickLabel',-5:5:30)
%     xlabel('Time [s]')
%     xlim([XLIM1 XLIM2])
    
%% Plot some variables
% figure
% subplot(3,3,1);
%     hold on
%     plot(nv.T, nv.out('K_e'), 'LineWidth', 1);
%     hold off
%     legend('W&C Neuron')     
%     ylabel('K_e');
%     xlim([XLIM1 XLIM2])
% subplot(3,3,2);
%     hold on
%     plot(nv.T, nv.out('R'), 'LineWidth', 1);
%     hold off
%     legend('W&C Neuron')
%     ylabel('R');
%     xlim([XLIM1 XLIM2])
% subplot(3,3,3);
%     hold on
%     plot(nv.T, nv.out('K_s'),'b');
%     hold off
%     legend('W&C Neuron')    
%     ylabel('K_s');
%     xlim([XLIM1 XLIM2])
% subplot(3,3,4);
%     hold on
%     plot(nv.T, nv.out('Na_sa'), 'LineWidth', 1);
%     hold off
%     legend('W&C Neuron')    
%     ylabel('Na_{sa}');
%     xlim([XLIM1 XLIM2])
% subplot(3,3,5);
%     hold on
%     plot(nv.T, nv.out('Na_d'), 'LineWidth', 1);
%     hold off
%     legend('W&C Neuron')    
%     ylabel('Na_d');
%     xlim([XLIM1 XLIM2])
% subplot(3,3,6);
%     hold on
%     plot(nv.T, nv.out('O2'), 'LineWidth', 1);
%     hold off
%     legend('W&C Neuron')    
%     ylabel('O_2');
%     xlim([XLIM1 XLIM2])
% subplot(3,3,7);
%     hold on
%     plot(nv.T, nv.out('NO_n'), 'LineWidth', 1);
%     hold off
%     legend('W&C Neuron')    
%     ylabel('NO_n');
%     xlim([XLIM1 XLIM2])
% subplot(3,3,8);
%     hold on
%     plot(nv.T, nv.out('Glu'), 'LineWidth', 1);
%     hold off
%     legend('W&C Neuron')
%     ylabel('Glu');
%     xlim([XLIM1 XLIM2])
% subplot(3,3,9);
%     hold on
%     plot(nv.T, nv.out('J_K_NEtoSC'),'b');
%     hold off
%     legend('W&C Neuron')     
%     ylabel('J_{KtoSC}');
%     xlim([XLIM1 XLIM2])
 
%% Plot variables compared with complex neuron model by Mathias et al.
% Load data from complex neuron model in NVU
% Mathias_Sq_nv = load('Dormanns_Squarpulse_withoutsecondpulse_nv.mat');
% Mathias_Sq_nv = struct2array(Mathias_Sq_nv);
% Ke16 = Mathias_Sq_nv.out('K_e');
% Ks16 = Mathias_Sq_nv.out('K_s');
% Nasa16 = Mathias_Sq_nv.out('Na_sa');
% Nad16 = Mathias_Sq_nv.out('Na_d');
% dKedt16 = Mathias_Sq_nv.out('J_K_NEtoSC')/11.5;
% O216 = Mathias_Sq_nv.out('O2');
% NO_n16 = Mathias_Sq_nv.out('NO_n');
% Glu16 = Mathias_Sq_nv.out('Glu');
% R16 = Mathias_Sq_nv.out('R');
% J_K16 = Mathias_Sq_nv.out('J_K_NEtoSC');
%  
% figure
% hold on
% plot(nv.T, nv.out('E_t'))
% hold on
% plot(nv.T, nv.out('E_t')-nv.out('I_t'))
% hold on
% plot(nv.T, nv.out('I_t'))
% hold off
% legend('E','E - I','I')
% ylabel('Fraction')
% xlim([XLIM1 XLIM2])
% figure
% plot(nv.T, nv.out('I_t'))
% ylabel('Proportion of inhibitory cells firing')
% xlim([XLIM1 XLIM2])
% 
% figure
% % subplot(3,3,1);
%     hold on
%     plot(nv.T, nv.out('K_e'), 'LineWidth', 1);
%     hold on
%     plot(nv.T, Ke16, 'LineWidth', 1);
%     hold off
%     legend('W&C Neuron','Elshin Neuron')     
%     ylabel('K_e');
%     xlim([XLIM1 XLIM2])
% % subplot(3,3,2);
% figure
%     hold on
%     plot(nv.T, nv.out('R'), 'LineWidth', 1);
%     hold on
%     plot(nv.T, R16)
%     hold off
%     legend('W&C Neuron','Elshin Neuron')
%     ylabel('R');
%     xlim([XLIM1 XLIM2])
% % subplot(3,3,3);
% figure
%     hold on
%     plot(nv.T, Ks16,'r')
%     hold on
%     plot(nv.T, nv.out('K_s'),'b');
%     hold off
%     legend('Elshin Neuron','W&C Neuron')    
%     ylabel('K_s');
%     xlim([XLIM1 XLIM2])
% % subplot(3,3,4);
% figure
%     hold on
%     plot(nv.T, nv.out('Na_sa'), 'LineWidth', 1);
%     hold on
%     plot(nv.T, Nasa16);
%     hold off
%     legend('W&C Neuron','Elshin Neuron')    
%     ylabel('Na_{sa}');
%     xlim([XLIM1 XLIM2])
% % subplot(3,3,5);
% figure
%     hold on
%     plot(nv.T, nv.out('Na_d'), 'LineWidth', 1);
%     hold on
%     plot(nv.T, Nad16);
%     hold off
%     legend('W&C Neuron','Elshin Neuron')    
%     ylabel('Na_d');
%     xlim([XLIM1 XLIM2])
% % subplot(3,3,6);
% figure
%     hold on
%     plot(nv.T, nv.out('O2'), 'LineWidth', 1);
%     hold on
%     plot(nv.T, O216);
%     hold off
%     legend('W&C Neuron','Elshin Neuron')    
%     ylabel('O_2');
%     xlim([XLIM1 XLIM2])
% % subplot(3,3,7);
% figure
%     hold on
%     plot(nv.T, nv.out('NO_n'), 'LineWidth', 1);
%     hold on
%     plot(nv.T, NO_n16)
%     hold off
%     legend('W&C Neuron','Elshin Neuron')    
%     ylabel('NO_n');
%     xlim([XLIM1 XLIM2])
% % subplot(3,3,8);
% figure
%     hold on
%     plot(nv.T, nv.out('Glu'), 'LineWidth', 1);
%     hold on
%     plot(nv.T, Glu16)
%     hold off
%     legend('W&C Neuron','Elshin Neuron')
%     ylabel('Glu');
%     xlim([XLIM1 XLIM2])
% % subplot(3,3,9);
% figure
%     hold on
%     plot(nv.T, J_K16,'r')
%     hold on
%     plot(nv.T, nv.out('J_K_NEtoSC'),'b');
%     hold off
%     legend('Elshin Neuron','W&C Neuron')     
%     ylabel('J_{KtoSC}');
%     xlim([XLIM1 XLIM2])


%% Additional BOLD signal plots
% figure;
%     plot(nv.T, BOLD_N, 'k', 'LineWidth', 1)
%     ylabel('BOLD (%)')
%     xlabel('Time [s]')
%     xlim([XLIM1 XLIM2])
%     p1=patch([NEURONAL_START NEURONAL_END NEURONAL_END NEURONAL_START],[-0.6 -0.6 1.2 1.2],'k');
%     set(p1,'FaceAlpha',0.1,'EdgeColor', 'none');
%     p2=patch([secondpulse secondpulse+secondpulselength ...
%     secondpulse+secondpulselength secondpulse],[-0.6 -0.6 1.2 1.2],'k');
%     set(p2,'FaceAlpha',0.1,'EdgeColor', 'none');
    
%% Plot Hemoglobin
% figure;
%     plot(nv.T, HBO_N, ':k', nv.T, HBR_N, '--k', nv.T, HBT_N, '-k', 'LineWidth', 1)
%     xlabel ('Time [s]')
%     xlim([XLIM1 XLIM2])
%     legend('HbO','HbR','HbT')
%     p1=patch([NEURONAL_START NEURONAL_END NEURONAL_END NEURONAL_START],[0.9 0.9 1.2 1.2],'k');
%     set(p1,'FaceAlpha',0.1,'EdgeColor', 'none');
%     p2=patch([secondpulse secondpulse+secondpulselength ...
%     secondpulse+secondpulselength secondpulse],[0.9 0.9 1.2 1.2],'k');
%     set(p2,'FaceAlpha',0.1,'EdgeColor', 'none');


%% Plot CBF
% figure;
% hold all;
% plot(nv.T, nv.out('CBF')./CBF_0, 'k', 'LineWidth', 1);
% ylabel('CBF [-]');
% xlim([XLIM1 XLIM2])
% ylim([0.97 1.4])
% title('Normalised CBF');
% xlabel('time (s)');
% p1=patch([100 116 116 100],[0.97 0.97 1.4 1.4],'k');
% set(p1,'FaceAlpha',0.1,'EdgeColor', 'none');
% p2=patch([124 125 125 124],[0.97 0.97 1.4 1.4],'k');
% set(p2,'FaceAlpha',0.1,'EdgeColor', 'none');


timeEnd = datetime('now');
fprintf('End time is %s\n', char(timeEnd));
