% Plot difference in power spectrum of the cytosolic calcium of cell 1 and
% 2, for different coupling D

clear;   

odeopts = odeset('RelTol', 1e-06, 'AbsTol', 1e-06, 'MaxStep', 0.1, 'Vectorized', 1);
[cmap]=buildcmap('kwr');

% Important constants:
J_PLC_1 = 0.3;
J_PLC_2 = 0.4;

startpulse = 0;       % ON: 10, OFF: 600, OFF/ON/OFF: 200
lengthpulse = 600;        % ON: 590, OFF: any, OFF/ON/OFF: 200

startofD = 0;
endofD = 0.15;
space = 0.001;

[cmap]=buildcmap('kwr');

fs = 10;    %sampling frequency (per second)
nfft = 2^13;

DD = startofD:space:endofD;

Z1 = zeros(length(DD),nfft/2-1);
Z2 = zeros(length(DD),nfft/2-1);
freq = fs*(1:nfft/2-1)/nfft;
i=1;
freq_cell1 = zeros(1,length(DD));
freq_cell2 = zeros(1,length(DD));

for D=startofD:space:endofD
    
    % SMC coupling
    D_Ca_i = D;
    D_IP3_i = 0;   % Optional IP3 coupling
    D_v_i = 0;     % Optional membrane potential coupling
    
    % EC coupling - optional!!
    D_Ca_j = 0;
    D_IP3_j = 0;   % Optional IP3 coupling
    D_v_j = 0;     % Optional membrane potential coupling
    
    nv = NVU_coupled(Astrocyte_1('startpulse', startpulse,'lengthpulse',lengthpulse), ...
    Astrocyte_2('startpulse', startpulse,'lengthpulse',lengthpulse), ...
    WallMechanics_1(), WallMechanics_2(), ...
    SMCEC_1('J_PLC_1', J_PLC_1, 'D_Ca_i', D_Ca_i, 'D_IP3_i', D_IP3_i, ...
    'D_v_i', D_v_i, 'D_Ca_j', D_Ca_j, 'D_IP3_j', D_IP3_j, 'D_v_j', D_v_j), ...
    SMCEC_2('J_PLC_2', J_PLC_2, 'D_Ca_i', D_Ca_i, 'D_IP3_i', D_IP3_i, ...
    'D_v_i', D_v_i, 'D_Ca_j', D_Ca_j, 'D_IP3_j', D_IP3_j, 'D_v_j', D_v_j), ...
    'odeopts', odeopts);

    nv.simulate();
    starttime = floor(300/nv.T(length(nv.T))*length(nv.T));
    endtime = length(nv.T);
    
    Ca_i_1 = nv.out('Ca_i_1'); s_i_1 = nv.out('s_i_1');
    Ca_i_2 = nv.out('Ca_i_2'); s_i_2 = nv.out('s_i_2');

    fftv1 = fft(Ca_i_1,nfft);
    absfftv1 = abs(fftv1); power1 = absfftv1(2:nfft/2);

    fftv2 = fft(Ca_i_2,nfft); 
    absfftv2 = abs(fftv2); power2 = absfftv2(2:nfft/2);

    Z1(i,:) = power1;
    Z2(i,:) = power2;
    
    i=i+1
end

figure;
imagesc(freq,(DD),Z1-Z2);
set(gca,'YDir','normal');
%xlim([0 0.8]);
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
