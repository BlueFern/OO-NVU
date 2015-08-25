% Power Spectrum for multiple D

clear;

odeopts = odeset('RelTol', 1e-04, 'AbsTol', 1e-04, 'MaxStep', 0.5, 'Vectorized', 1);

fs = 10;    %sampling frequency (per second)
nfft = 2^13;
endtime = 2000;
starttime = 400;

% Time
T = 0:1/fs:endtime;

% Important constants:
J_PLC_1 = 0.53;
J_PLC_2 = 0.58;

K_p = 9200;

startofD = 0;
endofD = 0.3;
space = 0.001;

[cmap]=buildcmap('kwr');

DD = startofD:space:endofD;
Z1 = zeros(length(DD),nfft/2-1);
Z2 = zeros(length(DD),nfft/2-1);
freq = fs*(1:nfft/2-1)/nfft;
i=1;
time_per_iteration = zeros(1,length(DD));
counter = 0;

num_iterations = length(DD);

fprintf('%d total iterations \n', num_iterations);

for D=startofD:space:endofD
    t = T;
    % SMC coupling
    D_Ca_i = D;
    D_IP3_i = 0;        % Optional IP3 coupling
    D_v_i = 0;          % Optional membrane potential coupling

    % EC coupling - optional!!
    D_Ca_j = 0;
    D_IP3_j = D_Ca_j;   % Optional IP3 coupling
    D_v_j = D_Ca_j;     % Optional membrane potential coupling
    
    nv = NVU_coupled_red( WallMechanics_1(), WallMechanics_2(), ...
    SMCEC_1_red('J_PLC_1', J_PLC_1, 'D_Ca_i', D_Ca_i, 'D_IP3_i', D_IP3_i, ...
    'D_v_i', D_v_i, 'D_Ca_j', D_Ca_j, 'D_IP3_j', D_IP3_j, 'D_v_j', D_v_j, 'K_p_1', K_p), ...
    SMCEC_2_red('J_PLC_2', J_PLC_2, 'D_Ca_i', D_Ca_i, 'D_IP3_i', D_IP3_i, ...
    'D_v_i', D_v_i, 'D_Ca_j', D_Ca_j, 'D_IP3_j', D_IP3_j, 'D_v_j', D_v_j, 'K_p_2', K_p), ...
    'odeopts', odeopts, 'T', t);

    nv.simulate();
    
    z1 = nv.out('Ca_i_1');
    z2 = nv.out('Ca_i_2');

    z1 = z1(fs*starttime:length(z1));
    z2 = z2(fs*starttime:length(z2));
    t = t(fs*starttime:length(t));

%     [maxcell1,indexmaxcell1] = max(abs(fft(z1-mean(z1),nfft)));
%     [maxcell2,indexmaxcell2] = max(abs(fft(z2-mean(z2),nfft)));
% 
%     freq_cell1(i) = freq(indexmaxcell1);
%     freq_cell2(i) = freq(indexmaxcell2);
    
    fftv1 = fft(z1,nfft);
    absfftv1 = abs(fftv1); power1 = absfftv1(2:nfft/2);

    fftv2 = fft(z2,nfft); 
    absfftv2 = abs(fftv2); power2 = absfftv2(2:nfft/2);

    Z1(i,:) = power1;
    Z2(i,:) = power2;
    
    i = i + 1;
    counter = counter + 1;
    amount_done = counter/num_iterations;
    percentage_done = amount_done*100;
    time_per_iteration(i) = nv.elapsedtime;
    time_left = mean(nonzeros(time_per_iteration))*(num_iterations-counter);
    min = floor(time_left/60);
    sec = rem(time_left, 60);
    
    fprintf('%.1f%% done, %d iterations left, ETA: %.f min %.f sec \n', percentage_done, num_iterations-counter, min, sec);
end

figure;
imagesc(freq,(DD),Z1-Z2);
set(gca,'YDir','normal');
xlim([0 0.5]);
ylim([startofD endofD]);
xlabel('Frequency (Hz)'); ylabel('D');
title('Difference Between the Cells');
colormap(cmap);
colorbar

% figure(101);
% plot(freq_cell1,DD,freq_cell2,DD,'LineWidth',1);
% ylabel('D');
% ylim([startofD endofD]);
% xlabel('Frequency');
% title('Frequency of oscillations');
% legend('Cell 1','Cell 2','Location','NorthEastOutside');

save powerspec9.mat
