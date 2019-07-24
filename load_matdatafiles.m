% Load experimental data from files and plot 

% clear all
% load './mat files/Hb_16_sec_Data.mat'
% figure(10)
% hold all
% plot(Time12, HbOMicromolarchange, 'r')
% plot(Time13, HbRMicromolarchange, 'b')
% plot(Time14, HbTMicromolarchange, 'g')
% legend('HbO', 'HbR', 'HbT')

clear all
load './mat files/Hb_whisker_pre7NI.mat'
figure(10)
hold all
plot(Time6-5, HbOPreMicromolarchange/100+1, '--r')
plot(Time7-5, HbRPreMicromolarchange/100+1, '--b')
plot(Time8-5, HbTPremicromolarchange/100+1, '--g')
legend('HbO', 'HbR', 'HbT')
title('Whisker response - Pre 7NI')
xlabel('Time [s]')
xlim([-2 20])

clear all
load './mat files/Hb_whisker_post7NI.mat'
figure(11)
hold all
plot(TIme9-5, HbOPost7NIMicromolarchange/100+1, '--r')
plot(Time10-5, HbRPost7NIMicromolarchange/100+1, '--b')
plot(Time11-5, HbTPost7NIMicromolarchange/100+1, '--g')
legend('HbO', 'HbR', 'HbT')
title('Whisker response - Post 7NI')
xlabel('Time [s]')
xlim([-2 20])

clear all
load './mat files/Hb_opto_pre7NI.mat'
figure(12)
hold all
plot(TimeHbO-5, HbO/100+1, '--r')
plot(TimeHbR-5, HbR/100+1, '--b')
plot(TimeHbT-5, HbT/100+1, '--g')
legend('HbO', 'HbR', 'HbT')
title('Opto response - Pre 7NI')
xlabel('Time [s]')
xlim([-2 20])

clear all
load './mat files/Hb_opto_post7NI.mat'
figure(13)
hold all
plot(TimeHbO-5, HbOMicromolarChange/100+1, '--r')
plot(TimeHbR-5, HbRMicromolarChange/100+1, '--b')
plot(TimeHbT-5, HbTMicromolarChange/100+1, '--g')
legend('HbO', 'HbR', 'HbT')
title('Opto response - Post 7NI')
xlabel('Time [s]')
xlim([-2 20])

%% L-NAME

clear all
load './mat files/Hb_whisker_preLNAME.mat'
figure(14)
hold all
plot(TimeHbO-5, HbO/100+1, '-.r')
plot(TimeHbR-5, HbR/100+1, '-.b')
plot(TimeHbT-5, HbT/100+1, '-.g')
legend('HbO', 'HbR', 'HbT')
title('Whisker response - Pre LNAME')
xlabel('Time [s]')
xlim([-2 20])

clear all
load './mat files/Hb_whisker_postLNAME.mat'
figure(15)
hold all
plot(Time3-5, HbOPostMicromolarchange/100+1, '-.r')
plot(Time4-5, HbRPostMicromolarchange/100+1, '-.b')
plot(Time5-5, HbTPostMicromolarchange/100+1, '-.g')
legend('HbO', 'HbR', 'HbT')
title('Whisker response - Post L-NAME')
xlabel('Time [s]')
xlim([-2 20])

clear all
load './mat files/Hb_opto_preLNAME.mat'
figure(16)
hold all
plot(Timeopt1-5, HbOMicromolarchange/100+1, '-.r')
plot(Timeopt2-5, HbRMicromolarchange/100+1, '-.b')
plot(TImeopt3-5, HbTMicromolarchange/100+1, '-.g')
legend('HbO', 'HbR', 'HbT')
title('Opto response - Pre L-NAME')
xlabel('Time [s]')
xlim([-2 20])

clear all
load './mat files/Hb_opto_postLNAME.mat'
figure(17)
hold all
plot(Timeopt4-5, HbOMicromolarchange/100+1, '-.r')
plot(Timeopt5-5, HbRMicromolarchange/100+1, '-.b')
plot(Timeopt6-5, HbTMicromolarchange/100+1, '-.g')
legend('HbO', 'HbR', 'HbT')
title('Opto response - Post L-NAME')
xlabel('Time [s]')
xlim([-2 20])

clear all

