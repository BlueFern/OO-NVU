%% Run script for NVU_DE_system - version 3

% To do:
% sort out other whisker cases
% create plot_all_variables.m

clear all;
timeStart = datetime('now');
fprintf('Start time is %s\n', char(timeStart));
tStart = tic;

% Add the path to the functions subfolder so they can be used (P, Q, Gamma functions)
addpath('./functions');

% Load indices, parameters and initial conditions
idx = all_indices();
p = all_parameters();
u0 = initial_conditions(idx);

% Set timespan and timestep for simulation [ms]
tspan = 0:p.dt:p.Tend;
numTimeSteps = length(tspan);

%% Load the experimental neural input data if needed
% Different Whiskerpad cases are described in all_parameters()
    
load Zheng2010_data.mat

if p.Whiskerpad == 1 || p.Whiskerpad == 2 || p.Whiskerpad == 4
    
    sum_neural_wh = zeros(size(neural_tim_vector)); % Initialise array
    for animal = 1:11
        for experiment = 1:10
            sum_neural_wh = sum_neural_wh + neural_data(:, p.ISI, p.stim, experiment, animal)'; % Sum all animal/trial data (110 total)
        end
    end
    mean_neural_wh = sum_neural_wh./110;  % Average the neural data

    neural_tim_vector_shifted = neural_tim_vector*1000 + p.startpulse - 20;    % Shift so stimulation begins at p.startpulse, -20 so initial spike isn't chopped off

    % Interpolate so there is data for all timesteps for NVU
    interp_neural_wh = interp1(neural_tim_vector_shifted, mean_neural_wh, tspan);
    interp_neural_wh(isnan(interp_neural_wh))=0;   % Remove NaNs     

    if p.double_pulse == 0  % Remove the second pulse if wanted
        t_idx = find(ismember(tspan, p.startpulse+p.lengthpulse+1e3));
        interp_neural_wh(t_idx:end) = 0;
    end
  
    % Locus coeruleus pathway - increases input stimulus after some period of time
    if p.LCpathway == 1
        time_1 = linspace(0, p.lengthpulse, 10000); 
        I_LC = 0.0004 * (time_1/1e3).^2 + 1;            % Chosen for the shape! Can be modified if needed
        time_1 = time_1 + p.startpulse;                 % Shift so starts at p.startpulse
        I_LC = interp1(time_1, I_LC, tspan);            % Make same size as tspan
        I_LC(isnan(I_LC)) = 1;                          % Replace NaNs with 1
    else
         I_LC = 1;
    end
    
    input_data = interp_neural_wh .* I_LC;   % Save input profile for use in P(t), Q(t)
        
    if p.Whiskerpad == 2 % Model with whiskerpad (Pinto et al) input - triangles for thalamic input T(t), overwrites input_data
        % Parameters
        Hz = 5;                   % [Hz or s^-1]                              % [individualpulses s^-1]
        peak_time = 1;            % [ms]                                      % Thalamic input increase time
        
        % Expressions
        msperpulse = 1000 ./ Hz;                                            % [ms individualpulse^-1]
        index_peak_time = peak_time + 1;                                    % Individual pulse increase length
        down_time = msperpulse - index_peak_time; %msperpulse - index_peak_time;                           % Individual pulse decrease length [ms], if we want the triangles to last only 15ms then use down_time = 15 - index_peak_time;
        rest_time = msperpulse - (down_time + index_peak_time);             % Rest time individual pulse
        numberOfPeriods = p.lengthpulse ./ msperpulse;                      % Number of individual pulses in neuronal stim

        % Make start and end time resting values
        start_time = p.startpulse;
        end_time = p.Tend - (p.startpulse + p.lengthpulse); %p.Tend+1;
        base_input = 0.0;           % [spikes/ms]                               % Base thalamic input
        startSignal = linspace(base_input, base_input, start_time/p.dt);
        endSignal = linspace(base_input, base_input, end_time/p.dt);

        % Construct one individual pulse cycle     
        firstamp = 1.07;             % [spikes/ms]                               % Thalamic input amplitude
        firstRiseSignal = linspace(base_input, firstamp, index_peak_time);
        firstFallingSignal = linspace(firstamp, base_input, down_time);
        restSignal = linspace(base_input, base_input, rest_time);

        firstCycle = [firstRiseSignal, firstFallingSignal(2:end-1), restSignal];

        % Now replicate this cycle several (numberOfPeriods) times.
        waveform = repmat(firstCycle, [1 numberOfPeriods]);

        % Interpolate so that it matches with the timestep p.dt 
        waveform_interp = interp1(linspace(0,p.lengthpulse,size(waveform,2)), waveform, 0:p.dt:p.lengthpulse); 
        
        % Add the bits before and after the stimulus
        total_waveform = [startSignal, waveform_interp, endSignal];
        
        % Multiply the waveforms by the Zheng input data (and LC pathway if used)
        input_data = total_waveform .* interp_neural_wh .* I_LC; 
        
    end
else
    input_data = 0; % Input data not used
end

clear interp_neural_wh mean_neural_wh neural_tim_vector neural_tim_vector_shifted sum_neural_wh neural_data experiment animal info

%% Run the model
options = odeset('RelTol', 1e-5, 'AbsTol', 1e-8, 'MaxStep', 1e3);
f = @(t,u) NVU_DE_system(t, u, idx, p, input_data);
[t,u] = ode15s(f,tspan,u0,options);

% Set the variable names for easy plotting and for use in all_algebraic_variables()
set_variable_names()
% Load all the fluxes for easy plotting
all_algebraic_variables()

% Output time to command window
tEnd = toc(tStart);
fprintf('Elapsed time is %d minutes and %f seconds\n',floor(tEnd/60),rem(tEnd,60));
timeEnd = datetime('now');
fprintf('End time is %s\n', char(timeEnd));

%% Calculate normalised CBF, CBV, HbR, CMRO2, HbT, HbO, %BOLD
preNeuronalStimTime = floor((p.startpulse)*numTimeSteps/p.Tend);
CBF_0 = CBF(preNeuronalStimTime);
CBV_0 = CBV(preNeuronalStimTime);
HBR_0 = HbR(preNeuronalStimTime);
CMRO2_0 = CMRO2(preNeuronalStimTime);

%Normalised variables
CBF_N = CBF./CBF_0;
CBV_N = CBV./CBV_0;
HBR_N = HbR./HBR_0;
CMRO2_N = CMRO2./CMRO2_0;

HBT_N = CBF_N .* HBR_N ./ CMRO2_N;                                      % Total hemoglobin (normalised), equivalent to HBR_N .* p.E_0 ./ OEF
HBO_N = (HBT_N - 1) - (HBR_N - 1) + 1;                                  % Oxyhemoglobin (normalised)
BOLD_N = 100 * p.V_0 * ( p.a_1 * (1 - HBR_N) - p.a_2 * (1 - CBV_N) );   % BOLD (percentage increase from 0)

%% Plot figures (with time in sec)
t_sec = t/1e3;                  % Convert time from ms to s
XLIM1 = p.startpulse/1e3 - 5;   % start plotting 5 sec before stimulus begins
XLIM2 = p.Tend/1e3;               % stop plotting at end of simulation time

% Plot the CBF from the Zheng experiments
if p.Whiskerpad == 1 || p.Whiskerpad == 2 || p.Whiskerpad == 3 || p.Whiskerpad == 4
    sum_cbf = zeros(size(cbf_tim_vector));
    for animal = 1:11
        for experiment = 1:10
            sum_cbf = sum_cbf+cbf_data(:,p.ISI,p.stim,experiment,animal)';
        end
    end
    mean_cbf = (sum_cbf./110) - 1;
    cbf_tim_vector_shifted = cbf_tim_vector * 1000 + p.startpulse; % Shift so stimulation begins at NEURONAL_START
    figure(300);
    hold all
    plot(cbf_tim_vector_shifted, mean_cbf, ':k', t_sec, (CBF - CBF_0) ./ CBF_0, '--k', 'LineWidth', 1);
    ylabel('\Delta CBF')
    xlabel('Time [s]')
    xlim([XLIM1 XLIM2])
    legend('Zheng experiment','Wilson & Cowan')
    title(['CBF with initial duration ' num2str(p.lengthpulse/1e3) 's, ISI ' num2str(p.actual_ISI) 's'] );
    p1=patch([p.startpulse p.startpulse+p.lengthpulse p.startpulse+p.lengthpulse p.startpulse]/1e3,[-0.05 -0.05 0.5 0.5],'k');
    set(p1,'FaceAlpha',0.1,'EdgeColor', 'none');
end
clear animal experiment sum_cbf

% Plot the 16 sec hemodynamics from the Berwick data
if p.lengthpulse == 16e3 && p.double_pulse == 0
    load scaledBerwickData.mat
    scaled_time_oxy = p.startpulse + scaled_time_oxy * 1e3;
    scaled_time_deoxy = p.startpulse + scaled_time_deoxy * 1e3;
    scaled_time_hbt = p.startpulse + scaled_time_hbt * 1e3;
    figure(301);
    hold all
    plot(scaled_time_oxy*1e-3, scaled_oxy, 'r', 'LineWidth', 1)
    plot(scaled_time_deoxy*1e-3, scaled_deoxy, 'b', 'LineWidth', 1)
    plot(scaled_time_hbt*1e-3, scaled_totaloxy, 'g', 'LineWidth', 1)
    xlabel('Time [s]')
    p1=patch([p.startpulse p.startpulse+p.lengthpulse p.startpulse+p.lengthpulse p.startpulse]/1e3,[0.97 0.97 1.05 1.05],'k');
    ylim([0.97 1.05]);
    set(p1,'FaceAlpha',0.05,'EdgeColor', 'none');
    title('Berwick experiment 16 sec: hemoglobin')
    xlim([XLIM1 XLIM2])
    legend('HbO','HbR','HbT');
end

% Plot the model hemodynamics
figure;
hold all;
plot(t_sec, HBO_N, 'r', 'LineWidth', 1)
plot(t_sec, HBR_N, 'b', 'LineWidth', 1)
plot(t_sec, HBT_N, 'g', 'LineWidth', 1)
xlim([XLIM1 XLIM2])
% ylim([0.87 1.2]);
xlabel('Time [s]')
% p1=patch([p.startpulse*1e-3 (p.startpulse + p.lengthpulse)*1e-3 (p.startpulse + p.lengthpulse)*1e-3 p.startpulse*1e-3],[0.7 0.7 1.3 1.3],'k');
% set(p1,'FaceAlpha',0.05,'EdgeColor', 'none');
legend('HbO','HbR','HbT');

% Plot the 20-HETE variables
% figure(23);
% subplot(3,2,1)
% hold all;
%     plot(t_sec, Ca_k);
%     xlabel('Time [s]'); ylabel('Ca_k [\muM]');
%     xlim([XLIM1 XLIM2])
% subplot(3,2,2)
% hold all;
%     plot(t_sec, H_i);
%     xlabel('Time [s]'); ylabel('H_i [\muM]');    
% %     ylim([0.12 0.24])
% xlim([XLIM1 XLIM2])
% subplot(3,2,3)
% hold all;
%     plot(t_sec, AA_k);
%     xlabel('Time [s]'); ylabel('AA_k [\muM]');
%     xlim([XLIM1 XLIM2])
% subplot(3,2,4)
% hold all;
%     plot(t_sec, AA_i);
%     xlabel('Time [s]'); ylabel('AA_i [\muM]');
%     xlim([XLIM1 XLIM2])
% subplot(3,2,5)
% hold all;
%     plot(t_sec, w_i);
%     xlabel('Time [s]'); ylabel('w_i [-]');
%     xlim([XLIM1 XLIM2])
% subplot(3,2,6)
% hold all;
%     plot(t_sec, R);
%     xlabel('Time [s]'); ylabel('R [\mum]');
%     xlim([XLIM1 XLIM2])
% subplot(4,2,7)
% hold all;
%     plot(t_sec, CMRO2);
%     xlabel('Time [s]'); ylabel('CMRO2');
%     xlim([XLIM1 XLIM2])

%% Plot anything else...
% 
% figure(3);
% % subplot(2,1,1)
% hold all
% plot(t_sec, (abs(E_t - I_t) - p.EImin)./ (p.EI_relative - p.EImin))
% ylabel('|E(t) - I(t)|')
% xlim([XLIM1 XLIM2])
% subplot(2,1,2)
% hold all
% if p.Whiskerpad == 2
%     plot(t_sec, T)
%     ylabel('T(t)')
% else
%     plot(t_sec, P)
%     ylabel('P(t)')    
% end
% xlim([XLIM1 XLIM2])


% figure(10);
% hold all
% plot(t_sec, (abs(E_t - I_t) - p.EImin)./ (p.EI_relative - p.EImin))
% ylabel('|E(t) - I(t)| - EI_{min} / EI_{max}')
% xlim([XLIM1 XLIM2])


% figure(12);
% hold all
% plot(t_sec, K_e)
% ylabel('K_e')
% xlim([XLIM1 XLIM2])

% figure(30);
% subplot(2,2,3)
% hold all
% plot(t_sec, OEF)
% ylabel('OEF')
% xlim([XLIM1 XLIM2])
% subplot(2,2,1)
% hold all
% plot(t_sec, f_in)
% ylabel('f_{in}')
% xlim([XLIM1 XLIM2])
% subplot(2,2,2)
% hold all
% plot(t_sec, CBF)
% ylabel('CBF')
% xlim([XLIM1 XLIM2])
% subplot(2,2,4)
% hold all
% plot(t_sec, CMRO2_N)
% ylabel('CMRO2')
% xlim([XLIM1 XLIM2])

figure(30)
subplot(3,3,1)
hold all
plot(t_sec, J_GABA_k)
ylabel('J_{GABA,k} [\muM/s]')
xlim([XLIM1 XLIM2])
subplot(3,3,2)
hold all
plot(t_sec, J_GABA_i)
ylabel('J_{GABA,i} [\muM/s]')
xlim([XLIM1 XLIM2])
subplot(3,3,3)
hold all
plot(t_sec, J_VOCC_i)
ylabel('J_{VOCC,i} [\muM/s]')
xlim([XLIM1 XLIM2])
subplot(3,3,4)
hold all
plot(t_sec, v_k)
ylabel('v_k [mV]')
xlim([XLIM1 XLIM2])
subplot(3,3,5)
hold all
plot(t_sec, v_i)
ylabel('v_i [mV]')
xlim([XLIM1 XLIM2])
subplot(3,3,6)
hold all
plot(t_sec, Ca_i)
ylabel('Ca_i [\muM]')
xlim([XLIM1 XLIM2])
subplot(3,3,7)
hold all
plot(t_sec, NO_n)
ylabel('NO_n [\muM]')
xlim([XLIM1 XLIM2])
subplot(3,3,8)
hold all
plot(t_sec, H_i)
ylabel('20-HETE [\muM]')
xlim([XLIM1 XLIM2])
subplot(3,3,9)
hold all
plot(t_sec, (CBF - CBF_0) ./ CBF_0)
ylabel('\Delta CBF')
xlim([XLIM1 XLIM2])


