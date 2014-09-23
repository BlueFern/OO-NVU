DATA = csvread(csvfilename);
time = DATA(:,length(DATA(1,:))-2);
close all
figure(1);
plot(time,DATA(:,flu.w_inf),'Linewidth',1.5)
xlabel('time in s','FontSize',12.0)
ylabel('w_inf [-]','FontSize',12.0)
title('equilibrium distribution of openings for the BK channel','FontSize',12.0)
set(gca,'fontsize',12)

figure(2)
plot(time,state(:,ind.w_k),'r','Linewidth',1.5)
xlabel('time in s','FontSize',12.0)
ylabel('w_k [-]','FontSize',12.0)
title('open BK channel probability','FontSize',12.0)
set(gca,'fontsize',12)

figure(3)
plot(time, DATA(:,flu.v_k)*1000,'Linewidth',1.5)
xlabel('time in s','FontSize',12.0)
ylabel('membrane voltage [V]','FontSize',12.0)
title('Astrocyte membrane voltage [V]','FontSize',12.0)
set(gca,'fontsize',12)

figure(4)
plot(time, DATA(:,flu.ck),'Linewidth',1.5)
xlabel('time in s','FontSize',12.0)
ylabel('Calcium concentration [uM]','FontSize',12.0)
title('calcium concentration','FontSize',12.0)
set(gca,'fontsize',12)

figure(5)
plot(time, DATA(:,flu.vh_3),'Linewidth',1.5)
xlabel('time [s]','FontSize',12.0)
ylabel('EET [uM]','FontSize',12.0)
title('EET concentration astrocyte [uM]','FontSize',12.0)
set(gca,'fontsize',12)

figure(6)
plot(time,DATA(:,flu.J_BK_k),'Linewidth',1.5)
xlabel('Time [s]','FontSize',12.0)
ylabel('K+ flux [uMs-1]','FontSize',12.0)
title('K+ flux through BK channel','FontSize',12.0)
set(gca,'fontsize',12)

figure(7)
plot(time,DATA(:,77+ind.K_p))
title('Potassium concentration in perivascular space')
xlabel('Time [s]')
ylabel('K_p  [\muM]')