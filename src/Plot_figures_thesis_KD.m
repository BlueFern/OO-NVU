%% Clear everything
clear all;
close all;
clc;

%% Construct NVU 
% odeopts = odeset('RelTol', 1e-03, 'AbsTol', 1e-03, 'MaxStep', 1, 'Vectorized', 1);
nv = NVU(Neuron(), ...
    Astrocyte(), ...
    WallMechanics(), ...
    SMCEC('J_PLC', 0.18)); % , ...
 %   'odeopts', odeopts);  % (change solver options)
nv.T = linspace(0, 1000, 2000); % NO pathway
% nv.T = linspace(0, 300, 1000); % NVU

% Run simulation
nv.simulate()

%% Decouple AC calcium:
nv.astrocyte.params.v_7 = 0;
% what else?

%% Switch off TRPV4 channel:
nv.astrocyte.params.switchBK = 1;
nv.astrocyte.params.reverseBK = 0; % V
nv.astrocyte.params.G_BK_k = 4.3e3; % pS (later converted to mho m^-2)
nv.astrocyte.params.trpv_switch = 0;
nv.astrocyte.params.epshalf_k = 0.16;
nv.astrocyte.params.Ca_4 = 0.15; % uM
nv.astrocyte.params.v_7 = -15e-3; % V

%% Switch off NO pathway:

c_w_i
K_act_i
K_2 / K_5


%% 1. plot: NVC overview
% yaxis{1} = '[K^+]_{sc} (\mu M)';        var_name{1} = 'K_s';
% yaxis{2} = 'v_k (V)';                   var_name{2} = 'v_k';
% yaxis{3} = 'J_{BK,k} (\mu M s^{-1})';   var_name{3} = 'J_BK_k';
% yaxis{4} = '[K^+]_{p} (\mu M)';         var_name{4} = 'K_p';
% yaxis{5} = 'J_{KIR,i} (\mu M s^{-1})';  var_name{5} = 'J_KIR_i';
% yaxis{6} = 'v_i (mV)';                  var_name{6} = 'v_i';
% yaxis{7} = 'J_{VOCC,i} (\mu M s^{-1})'; var_name{7} = 'J_VOCC_i';
% yaxis{8} = '[Ca^{2+}]_i (\mu M)';       var_name{8} = 'Ca_i';
% yaxis{9} = 'F (dim.less)';              var_name{9} = 'F_r';
% yaxis{10} = 'r (m)';                    var_name{10} = 'R';
% close all

figure(99); 
set(gcf,'Name','NVC overview')
% set(gcf,'Position', [24 62 1616 904], ...
%         'PaperPosition', [0.634517 6.34517 20.3046 15.2284]);
set(gcf,'Position', [24 62 904 1616], ...
        'PaperPosition', [0.634517 6.34517 20.3046 15.2284]);
    
subplot(5,2,1)
plot(nv.T, 1e-3*nv.out('K_s'))
% xlabel('Time (s)')
ylabel('[K^+]_{sc} (mM)') % converted from \mu
% ylabel('\\centering[K$\^+$]$\_\{\\text\{sc\}\}$ (\\textmu M)')
% grid on; grid minor
hold all
   
subplot(5,2,2)
plot(nv.T, 1e3*nv.out('v_k'))
% xlabel('Time (s)')
ylabel('v_k (mV)')
% grid on; grid minor
hold all

% subplot(5,2,3)
% plot(nv.T, 1e-3*nv.out('J_BK_p')) % K+ flux into the PVS 
% % xlabel('Time (s)')
% ylabel('J_{BK,p} (mM s^{-1})') % converted from \mu
% grid on; grid minor
% hold all

subplot(5,2,3)
plot(nv.T, nv.out('J_BK_k')) % K+ flux into the PVS 
% xlabel('Time (s)')
ylabel('J_{BK,k} (\muM s^{-1})')
% grid on; grid minor
hold all

subplot(5,2,4)
plot(nv.T, 1e-3*nv.out('K_p'))
% xlabel('Time (s)')
ylabel('[K^+]_{p} (mM)')
% grid on; grid minor
hold all

subplot(5,2,5)
plot(nv.T, nv.out('J_KIR_i'))
% xlabel('Time (s)')
ylabel('J_{KIR,i} (\muM s^{-1})')
% grid on; grid minor
hold all

subplot(5,2,6)
plot(nv.T, nv.out('v_i'))
% xlabel('Time (s)')
ylabel('v_i (mV)')
% grid on; grid minor
hold all

subplot(5,2,7)
plot(nv.T, 1e3*nv.out('J_VOCC_i'))
% xlabel('Time (s)')
ylabel('J_{VOCC,i} (nM s^{-1})') % converted from \mu
% grid on; grid minor
hold all

subplot(5,2,8)
plot(nv.T, nv.out('Ca_i'))
% xlabel('Time (s)')
ylabel('[Ca^{2+}]_i (\muM)')
% grid on; grid minor
hold all

subplot(5,2,9)
plot(nv.T, nv.out('F_r'))
xlabel('Time (s)')
ylabel('F_r (dim.less)')
% grid on; grid minor
hold all

subplot(5,2,10)
plot(nv.T, 1e6*nv.out('R'))
xlabel('Time (s)')
ylabel('r (\mum)')
% grid on; grid minor
hold all
% matlab2tikz('~/Dropbox/1PhD_UCCanterbury/PhD_Thesis_KD/tikz_figures/all2.tikz', 'height', '\figureheight', 'width', '\figurewidth');

%%
figure(1); 
plot(nv.T, 1e-3*nv.out('K_s'))
% xlabel('Time (s)')
ylabel('[K^+]_{sc} (mM)') % converted from \mu
hold all
%matlab2tikz('~/Dropbox/1PhD_UCCanterbury/PhD_Thesis_KD/tikz_figures/a.tikz', 'height', '\figureheight', 'width', '\figurewidth');
   
figure(2)
plot(nv.T, 1e3*nv.out('v_k'))
% xlabel('Time (s)')
ylabel('v_k (mV)')
hold all
%matlab2tikz('~/Dropbox/1PhD_UCCanterbury/PhD_Thesis_KD/tikz_figures/b.tikz', 'height', '\figureheight', 'width', '\figurewidth');

figure(3)
plot(nv.T, nv.out('J_BK_k')) % K+ flux into the PVS 
% xlabel('Time (s)')
ylabel('J_{BK,k} (\muM s^{-1})')
hold all
%matlab2tikz('~/Dropbox/1PhD_UCCanterbury/PhD_Thesis_KD/tikz_figures/c.tikz', 'height', '\figureheight', 'width', '\figurewidth');

figure(4)
plot(nv.T, 1e-3*nv.out('K_p'))
axis([0 300 0 15])
% xlabel('Time (s)')
ylabel('[K^+]_p (mM)')
hold all
%matlab2tikz('~/Dropbox/1PhD_UCCanterbury/PhD_Thesis_KD/tikz_figures/d.tikz', 'height', '\figureheight', 'width', '\figurewidth');

figure(5)
plot(nv.T, nv.out('J_KIR_i'))
axis([0 300 0 8e-2])
% xlabel('Time (s)')
ylabel('J_{KIR,i} (\muM s^{-1})')
hold all
%matlab2tikz('~/Dropbox/1PhD_UCCanterbury/PhD_Thesis_KD/tikz_figures/e.tikz', 'height', '\figureheight', 'width', '\figurewidth');

figure(6)
plot(nv.T, nv.out('v_i'))
axis([0 300 -60 -20])
% xlabel('Time (s)')
ylabel('v_i (mV)')
hold all
%matlab2tikz('~/Dropbox/1PhD_UCCanterbury/PhD_Thesis_KD/tikz_figures/f.tikz', 'height', '\figureheight', 'width', '\figurewidth');

figure(7)
plot(nv.T, 1e3*nv.out('J_VOCC_i'))
% xlabel('Time (s)')
ylabel('J_{VOCC,i} (nM s^{-1})') % converted from \mu
hold all
%matlab2tikz('~/Dropbox/1PhD_UCCanterbury/PhD_Thesis_KD/tikz_figures/g.tikz', 'height', '\figureheight', 'width', '\figurewidth');

figure(8)
plot(nv.T, nv.out('Ca_i'))
axis([0 300 0 1])
% xlabel('Time (s)')
ylabel('[Ca^{2+}]_i (\muM)')
hold all
%matlab2tikz('~/Dropbox/1PhD_UCCanterbury/PhD_Thesis_KD/tikz_figures/h.tikz', 'height', '\figureheight', 'width', '\figurewidth');

figure(9)
plot(nv.T, nv.out('F_r'))
xlabel('Time (s)')
ylabel('F_r (dim.less)')
hold all
%matlab2tikz('~/Dropbox/1PhD_UCCanterbury/PhD_Thesis_KD/tikz_figures/i.tikz', 'height', '\figureheight', 'width', '\figurewidth');

figure(10)
plot(nv.T, 1e6*nv.out('R'))
axis([0 300 16 26])
xlabel('Time (s)')
ylabel('r (\mum)')
hold all
%matlab2tikz('~/Dropbox/1PhD_UCCanterbury/PhD_Thesis_KD/tikz_figures/j.tikz', 'height', '\figureheight', 'width', '\figurewidth');


% matlab2tikz('~/Dropbox/1PhD_UCCanterbury/PhD_Thesis_KD/tikz_figures/NVU_FIG1.tikz', 'height', '\figureheight', 'width', '\figurewidth');

    
% subplot(3, 2, 1)
% plot(nv.T, nv.out(''))
% xlabel('')
% ylabel('')
% 
% set(gcf,'Name','SMC fluxes')
% set(gcf,'Position', [24 62 1616 904],...
%         'PaperPosition', [0.634517 6.34517 20.3046 15.2284]...
%         );
% for j = [1:3 5:16]
% subplot(4,4,j)
% plot(time,DATA(:,smcoff+j))
% % h101(j+2) = gca();
% xlabel('time in s')
% ylabel(SMCtitle{j})
% title(SMCtitle{j})
% hold all
% end


%% 2. plot: Contribution of different Ca2+ fluxes towards SMC cytosol [Ca2+]



%% NO pathway plots
% nv.T = linspace(0, 1000, 2000);
% nv.neuron.params.startpulse = 200;
% nv.neuron.params.lengthpulse = 200;
% nv.astrocyte.params.startpulse = 200;
% nv.astrocyte.params.lengthpulse = 200;
% nv.simulate()

figure(7); 
set(gcf,'Name','NO pathway')
% set(gcf,'Position', [24 62 1616 904], ...
%         'PaperPosition', [0.634517 6.34517 20.3046 15.2284]);
set(gcf,'Position', [24 62 904 1616], ...
        'PaperPosition', [0.634517 6.34517 20.3046 15.2284]);
    
subplot(3,2,1)
plot(nv.T, nv.out('nNOS_act_n'),nv.T, nv.out('eNOS_act_j'))
% xlabel('Time (s)')
ylabel('[NOS]')
ylim([0 1.5])
% ylabel('\\centering[K$\^+$]$\_\{\\text\{sc\}\}$ (\\textmu M)')
grid on; grid minor
hold all

subplot(3,2,2)
plot(nv.T, nv.out('NO_n'),nv.T, nv.out('NO_k'),nv.T, nv.out('NO_i'),nv.T, nv.out('NO_j'))
% xlabel('Time (s)')
ylabel('NO')
% ylabel('\\centering[K$\^+$]$\_\{\\text\{sc\}\}$ (\\textmu M)')
grid on; grid minor
hold all

subplot(3,2,3)
plot(nv.T, nv.out('cGMP_i'))
% xlabel('Time (s)')
ylabel('cGMP')
% ylabel('\\centering[K$\^+$]$\_\{\\text\{sc\}\}$ (\\textmu M)')
grid on; grid minor
hold all

subplot(3,2,4)
plot(nv.T, nv.out('K_2'))
% xlabel('Time (s)')
ylabel('K2')
% ylabel('\\centering[K$\^+$]$\_\{\\text\{sc\}\}$ (\\textmu M)')
grid on; grid minor
hold all

subplot(3,2,5)
plot(nv.T, nv.out('w_i'))
% xlabel('Time (s)')
ylabel('wi')
% ylabel('\\centering[K$\^+$]$\_\{\\text\{sc\}\}$ (\\textmu M)')
grid on; grid minor
hold all

subplot(3,2,6)
plot(nv.T, nv.out('R'))
% xlabel('Time (s)')
ylabel('R')
% ylabel('\\centering[K$\^+$]$\_\{\\text\{sc\}\}$ (\\textmu M)')
grid on; grid minor
hold all

%% 
figure(11)
subplot(2,2,1)
plot(nv.T, nv.out('E_b'))

subplot(2,2,2)
plot(nv.T, nv.out('E_6c'))
grid minor;

subplot(2,2,3)
plot(nv.T, nv.out('E_5c'))
grid on;

%matlab2tikz('~/Dropbox/1PhD_UCCanterbury/PhD_Thesis_KD/tikz_figures/NVU_FIG_test.tikz', 'height', '\figureheight', 'width', '\figurewidth');

%% 
% figure(12)
% subplot(2,2,1)
% reduce_plot(nv.T, nv.out('E_b'))
% 
% subplot(2,2,2)
% reduce_plot(nv.T, nv.out('E_6c'))
% 
% subplot(2,2,3)
% reduce_plot(nv.T, nv.out('E_5c'))
% 
% matlab2tikz('~/Dropbox/1PhD_UCCanterbury/PhD_Thesis_KD/tikz_figures/NVU_FIG_test_reduced.tikz', 'height', '\figureheight', 'width', '\figurewidth');

%% J_PLC = 0.18 / 0.4
    close all
    nv.T = linspace(0, 300, 1000);
%     nv.neuron.params.startpulse = 100;
%     nv.neuron.params.lengthpulse = 100;
%     nv.astrocyte.params.startpulse = 100;
%     nv.astrocyte.params.lengthpulse = 100;
    
for i = [0.18 ,0.4]
    nv.smcec.params.J_PLC = i;
    nv.simulate()

    figure(1); 
    plot(nv.T, 1e-3*nv.out('K_s'))
    % xlabel('Time (s)')
    ylabel('[K^+]_{sc} (mM)') % converted from \mu
    hold all
    %matlab2tikz('~/Dropbox/1PhD_UCCanterbury/PhD_Thesis_KD/tikz_figures/a2.tikz', 'height', '\figureheight', 'width', '\figurewidth');

    figure(2)
    plot(nv.T, 1e3*nv.out('v_k'))
    % xlabel('Time (s)')
    ylabel('v_k (mV)')
    hold all
    %matlab2tikz('~/Dropbox/1PhD_UCCanterbury/PhD_Thesis_KD/tikz_figures/b2.tikz', 'height', '\figureheight', 'width', '\figurewidth');

    figure(3)
    plot(nv.T, nv.out('J_BK_k')) % K+ flux into the PVS 
    % xlabel('Time (s)')
    ylabel('J_{BK,k} (\muM s^{-1})')
    hold all
    %matlab2tikz('~/Dropbox/1PhD_UCCanterbury/PhD_Thesis_KD/tikz_figures/c2.tikz', 'height', '\figureheight', 'width', '\figurewidth');

    figure(4)
    plot(nv.T, 1e-3*nv.out('K_p'))
    axis([0 300 0 15])
    % xlabel('Time (s)')
    ylabel('[K^+]_p (mM)')
    hold all
    %matlab2tikz('~/Dropbox/1PhD_UCCanterbury/PhD_Thesis_KD/tikz_figures/d2.tikz', 'height', '\figureheight', 'width', '\figurewidth');

    figure(5)
    plot(nv.T, nv.out('J_KIR_i'))
    axis([0 300 0 8e-2])
    % xlabel('Time (s)')
    ylabel('J_{KIR,i} (\muM s^{-1})')
    hold all
    %matlab2tikz('~/Dropbox/1PhD_UCCanterbury/PhD_Thesis_KD/tikz_figures/e2.tikz', 'height', '\figureheight', 'width', '\figurewidth');

    figure(6)
    plot(nv.T, nv.out('v_i'))
    axis([0 300 -60 -20])
    % xlabel('Time (s)')
    ylabel('v_i (mV)')
    hold all
    %matlab2tikz('~/Dropbox/1PhD_UCCanterbury/PhD_Thesis_KD/tikz_figures/f2.tikz', 'height', '\figureheight', 'width', '\figurewidth');

    figure(7)
    plot(nv.T, 1e3*nv.out('J_VOCC_i'))
    % xlabel('Time (s)')
    ylabel('J_{VOCC,i} (nM s^{-1})') % converted from \mu
    hold all
    %matlab2tikz('~/Dropbox/1PhD_UCCanterbury/PhD_Thesis_KD/tikz_figures/g2.tikz', 'height', '\figureheight', 'width', '\figurewidth');

    figure(8)
    plot(nv.T, nv.out('Ca_i'))
    axis([0 300 0 1])
    % xlabel('Time (s)')
    ylabel('[Ca^{2+}]_i (\muM)')
    hold all
    %matlab2tikz('~/Dropbox/1PhD_UCCanterbury/PhD_Thesis_KD/tikz_figures/h2.tikz', 'height', '\figureheight', 'width', '\figurewidth');

    figure(9)
    plot(nv.T, nv.out('F_r'))
    xlabel('Time (s)')
    ylabel('F_r (dim.less)')
    hold all
    %matlab2tikz('~/Dropbox/1PhD_UCCanterbury/PhD_Thesis_KD/tikz_figures/i2.tikz', 'height', '\figureheight', 'width', '\figurewidth');

    figure(10)
    plot(nv.T, 1e6*nv.out('R'))
    axis([0 300 16 26])
    xlabel('Time (s)')
    ylabel('r (\mum)')
    hold all
    %matlab2tikz('~/Dropbox/1PhD_UCCanterbury/PhD_Thesis_KD/tikz_figures/j2.tikz', 'height', '\figureheight', 'width', '\figurewidth');
end


%% Scaling results HP1 and POWER7

no_processors = [32, 16, 8, 4, 2, 1];
weak_POWER7 = [12.2344722222 11.8798611111 11.7373611111 11.5920277778 11.3818055556 11.1873333333];
weak_HP1 = [5.2310833333 5.1254444444 5.0393888889 5.0174722222 4.9570555556 4.8980277778];
strong_POWER7 = [1.3452055556 1.8976027778 3.2069166667 5.8918888889 11.3818055556 22.1152222222];
strong_HP1 = [0.5459305556 0.8061222222 1.3058638889 2.2926611111 4.9570555556 8.1811111111];
% ideal_strong = [1/32, 1/16, 1/8, 1/4, 1/2, 1];


% figure; loglog(no_processors,weak_POWER7,':+',no_processors,weak_HP1,':o');
figure; semilogx(no_processors,weak_POWER7,':+',no_processors,weak_HP1,':o');
axis([1 100 0 15])
%matlab2tikz('~/Dropbox/1PhD_UCCanterbury/PhD_Thesis_KD/tikz_figures/weak_scaling.tikz', 'height', '\figureheight', 'width', '\figurewidth');


figure; loglog(no_processors,strong_POWER7,':+',no_processors,strong_HP1,':o');
%matlab2tikz('~/Dropbox/1PhD_UCCanterbury/PhD_Thesis_KD/tikz_figures/strong_scaling.tikz', 'height', '\figureheight', 'width', '\figurewidth');

% figure; loglog(no_processors,ideal_strong);
% axis([0 32])
