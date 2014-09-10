time = DATA(:,length(DATA(1,:))-2);
close all
figure(1);
plot(time,DATA(:,acoff+flu.w_inf),'Linewidth',1.5)
xlabel('time in s','FontSize',12.0)
ylabel('w_inf [-]','FontSize',12.0)
title('equilibrium distribution of openings of the BK channel','FontSize',12.0)
set(gca,'fontsize',12)

figure(2)
plot(time, state(:,ind.w_i),'r','Linewidth',1.5)
xlabel('time in s','FontSize',12.0)
ylabel('w_i [-]','FontSize',12.0)
title('open probability K^+ channel','FontSize',12.0)
set(gca,'fontsize',12)

figure(3)
plot(time, -state(:,ind.eetk),'Linewidth',1.5)
xlabel('time in s','FontSize',12.0)
ylabel('EET concentration [uM]','FontSize',12.0)
title('EET concentration','FontSize',12.0)
set(gca,'fontsize',12)

figure(4)
plot(time, -state(:,ind.ck),'Linewidth',1.5)
xlabel('time in s','FontSize',12.0)
ylabel('Calcium concentration [uM]','FontSize',12.0)
title('calcium concentration','FontSize',12.0)
set(gca,'fontsize',12)

figure(5)
plot(time, state(:,ind.N_K_k),'Linewidth',1.5)
xlabel('time in s','FontSize',12.0)
ylabel('Potassium concentration astrocyte [uM]','FontSize',12.0)
title('K+ [uM]','FontSize',12.0)
set(gca,'fontsize',12)

figure(6)
plot(time,DATA(:,acoff+flu.ck),'Linewidth',1.5)
xlabel('Time [s]','FontSize',12.0)
ylabel('v_k [mV]','FontSize',12.0)
title('Membrane Potential of the astrocyte','FontSize',12.0)
set(gca,'fontsize',12)