clear;

odeopts = odeset('RelTol', 1e-04, 'AbsTol', 1e-04, 'MaxStep', 0.5, 'Vectorized', 1);

fs = 10;    %sampling frequency (per second)
nfft = 2^13;
endtime = 2000;
starttime = 400;

% Time
t = 0:1/fs:endtime;

% Important constants:
J_PLC_1 = 0.25;
J_PLC_2 = 0.3;

% SMC coupling
D_Ca_i = 2;
D_IP3_i = 0;        % Optional IP3 coupling
D_v_i = 0;          % Optional membrane potential coupling

% EC coupling - optional!!
D_Ca_j = 0;
D_IP3_j = D_Ca_j;   % Optional IP3 coupling
D_v_j = D_Ca_j;     % Optional membrane potential coupling

K_p = 9200;

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

freq = fs*(2:nfft/2-1)/nfft;

% [maxcell1,indexmaxcell1] = max(abs(fft(z1-mean(z1),nfft)));
% [maxcell2,indexmaxcell2] = max(abs(fft(z2-mean(z2),nfft)));

% freqcell1 = freq(indexmaxcell1)
% freqcell2 = freq(indexmaxcell2)
% period1 = 1/freqcell1
% period2 = 1/freqcell2


fftz1 = fft(z1,nfft);
absfftz1 = abs(fftz1); 

fftz2 = fft(z2,nfft);
absfftz2 = abs(fftz2); 

figure(40);
plot(fs*(2:nfft/2-1)/nfft, absfftz1(3:nfft/2),'r','LineWidth',1); %Ignore first point - very high
xlim([0 1]);
title('Power Spectrum of the Cytosolic Calcium in Cell 1');
xlabel('Frequency (Hz)');
ylabel('Power');

figure(41);
plot(fs*(2:nfft/2-1)/nfft, absfftz2(3:nfft/2),'r','LineWidth',1); %Ignore first point - very high
xlim([0 1]);
title('Power Spectrum of the Cytosolic Calcium in Cell 2');
xlabel('Frequency (Hz)');
ylabel('Power');

figure(43);
plot(fs*(2:nfft/2-1)/nfft, absfftz1(3:nfft/2)-absfftz2(3:nfft/2));
xlim([0 1]);
title('Difference between the cells');
xlabel('Frequency (Hz)');
ylabel('Power');