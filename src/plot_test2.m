% close all
figure(8)
plot(time,(DATA(:,acoff+flu.J_BK_k)./DATA(:,stoff+ind.R_k)))
title('K^+ flux through the BK channel')
xlabel('Time [s]')
ylabel('K^+ flux [\muM/s]')

figure(9)
plot(time,DATA(:,stoff+ind.R_k))

figure(10)
plot(time,DATA(:,flu.J_K_k))


% hold on
% DATA = csvread(csvfilename);
% plot(t,state(:,flu.Na_k),'r','Linewidth',1.5)
% figure(2)
% plot(t,DATA(:,flu.ck),'Linewidth',1.5)
% 
% figure(3)
% plot(t,DATA(:,flu.K_k),'k', 'Linewidth',1.5)
% 
% xlabel('Time [s]','FontSize',12.0)
% ylabel('v_k [mV]','FontSize',12.0)
% title('Membrane Potential of the astrocyte','FontSize',12.0)
% set(gca,'fontsize',12)


