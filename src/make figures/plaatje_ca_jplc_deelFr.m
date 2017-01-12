%graph you can make after running figure_ca_jplc_deel2

lijst=[1 3 5];
legendlist=[];
for i=[1:length(lijst)]
    nummer=lijst(i);
    figure(11)
    
    kleur=['k','b','r','g'];
    hold on
    plot(simulatie.J_PLC_value(nummer,:),simulatie.Frmaximum(nummer,:), 'color',kleur(i),'LineWidth',2.5)
    plot(simulatie.J_PLC_value(nummer,:),simulatie.Frminimum(nummer,:), 'color',kleur(i),'LineWidth',2.5)
 %   plot(simulatie.J_PLC_value(nummer,:),simulatie.NOmean(nummer,:), 'color',kleur(i),'LineWidth',2.5)
 %   myfit=polyfit(simulatie.J_PLC_value(nummer,:),simulatie.radiusmean(nummer,:)*0.1*1e6,4)
    
    X22(1,:)=0.2:0.25/6:0.45;
    X22(2,:)=0.25:0.275/6:0.525;
    X22(3,:)=0.1:0.185/6:0.285;
%    Y2=polyval(myfit,X22(i,:));
   % plot(X22(i,:),Y2, '--', 'color',kleur(i),'LineWidth',2.5)
    xlabel('J PLC [\muM/s]', 'FontSize', 26); 
    ylabel('Fraction of attached myosin [-]', 'FontSize', 26)
    box on
    %title('Membrame potential with different pressure')    
    set(gca,'FontSize',30,'LineWidth',2)
    
    tix=get(gca,'ytick')';
   % set(gca,'yticklabel',num2str(tix,'%.2f'))
    simulatie.description(nummer,:);
    if length(legendlist)==0;
        legendlist=[simulatie.description(nummer,:)];
    else
        legendlist=[legendlist ',' simulatie.description(nummer,:)];
    end
    
end

legenda = strsplit(legendlist,',');
legend(legenda);