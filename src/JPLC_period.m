% Plot period of calcium with JPLC as parameter

clear;   

odeopts = odeset('RelTol', 1e-06, 'AbsTol', 1e-06, 'MaxStep', 0.1, 'Vectorized', 1);

% Signal on:
startofjplc = 0.27;
endofjplc = 0.56;

%Signal off:
% startofjplc = 0.23;
% endofjplc = 0.5;

space = 0.005;

i=1;
fs = 10;    %sampling frequency (per second)
nfft = 2^13;
jplc = startofjplc:space:endofjplc;
Z1 = zeros(length(jplc),nfft/2-1);
freq = fs*(1:nfft/2-1)/nfft;
freq_cell1 = zeros(1,length(jplc));

for JPLC=startofjplc:space:endofjplc
    
    % Signal on
    nv = NVU(Astrocyte('startpulse', 10, 'lengthpulse', 590), ...
    WallMechanics(), ...
    SMCEC('J_PLC', JPLC), ...
    'odeopts', odeopts);
     
    % Signal off
%     nv = NVU(Astrocyte('startpulse', 600), ...
%     WallMechanics(), ...
%     SMCEC('J_PLC', JPLC), ...
%     'odeopts', odeopts);



    nv.simulate();
    
    starttime = floor(200/nv.T(length(nv.T))*length(nv.T));
    endtime = length(nv.T);
    
    z1 = nv.out('Ca_i');
    z1 = z1(starttime:endtime);
    
    [maxcell1,indexmaxcell1] = max(abs(fft(z1-mean(z1),nfft)));

    freq_cell1(i) = freq(indexmaxcell1);
    
%     fftv1 = fft(Ca_i_1,nfft);
%     absfftv1 = abs(fftv1); power1 = absfftv1(2:nfft/2);
% 
%     fftv2 = fft(Ca_i_2,nfft); 
%     absfftv2 = abs(fftv2); power2 = absfftv2(2:nfft/2);

%     Z1(i,:) = power1;
%     Z2(i,:) = power2;

    i = i+1
end

figure(102);
hold on
%plot(jplc,1./freq_cell1,'r','LineWidth',1); % Signal on
plot(jplc,1./freq_cell1,'b','LineWidth',1); % Signal off
hold off
ylabel('Period');
%xlim([0.2 0.6]);
%ylim([0 40]);
xlim([0.22 0.29]);
ylim([14 18]);
xlabel('J_{PLC}');
title('Period of oscillations');