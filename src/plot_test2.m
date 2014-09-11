close all
figure(1)
hold on
plot(time,state(:,acoff+flu.Na_k),'r','Linewidth',1.5)
figure(2)
plot(time,DATA(:,acoff+flu.ck),'Linewidth',1.5)

figure(3)
plot(time,DATA(:,acoff+flu.K_k),'k', 'Linewidth',1.5)

xlabel('Time [s]','FontSize',12.0)
ylabel('v_k [mV]','FontSize',12.0)
title('Membrane Potential of the astrocyte','FontSize',12.0)
set(gca,'fontsize',12)


a=DATA(:,acoff+flu.ck)
b=DATA(:,acoff+flu.K_k)
