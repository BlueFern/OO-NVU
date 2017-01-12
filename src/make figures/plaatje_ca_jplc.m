%some graphs you can make after running figure_ca_jplc_deel2

lijst=[1 3 5];
legendlist=[];
for i=[1:length(lijst)]
    nummer=lijst(i);
    figure(1)
    
    kleur=['k','b','r','g'];
    hold on
    plot(simulatie.J_PLC_value(nummer,:),simulatie.maximumlist(nummer,:), 'color',kleur(i),'LineWidth',2.5)
    plot(simulatie.J_PLC_value(nummer,:),simulatie.minimumlist(nummer,:), 'color',kleur(i),'LineWidth',2.5)
    xlabel('J PLC [\muM/s]', 'FontSize', 22); 
    ylabel('Ca2+ concentration SMC [\muM]', 'FontSize', 22)
    box on
    title('Position Hopf bifurcations with different pressure')
    set(gca,'FontSize',22,'LineWidth',2)
    
    simulatie.description(nummer,:);
    if length(legendlist)==0;
        legendlist=[simulatie.description(nummer,:)];
    else
        legendlist=[legendlist ',' simulatie.description(nummer,:)];
    end
    
end

legenda = strsplit(legendlist,',');
legend(legenda);

%lijst=[2];
legendlist=[];
for i=[1:length(lijst)]
    nummer=lijst(i);
    figure(2)
    subplot(2,1,1)
    kleur=['k','b','r','g'];
    hold on
    plot(simulatie.J_PLC_value(nummer,:),simulatie.maximumlist(nummer,:).^4 ./ (simulatie.maximumlist(nummer,:).^4 + 0.9^4), 'color',kleur(i))
    plot(simulatie.J_PLC_value(nummer,:),simulatie.minimumlist(nummer,:).^4 ./ (simulatie.minimumlist(nummer,:).^4 + 0.9^4), 'color',kleur(i))
    xlabel('J PLC [\muM/s]'); title('calciumpart [-]')
    simulatie.description(nummer,:);
    if length(legendlist)==0;
        legendlist=[simulatie.description(nummer,:)];
    else
        legendlist=[legendlist ',' simulatie.description(nummer,:)];
    end
    
end
legenda = strsplit(legendlist,',');
legend(legenda);

%lijst=[2];
legendlist=[];
for i=[1:length(lijst)]
    nummer=lijst(i);
    figure(2)
    subplot(2,1,2)
    kleur=['k','b','r','g'];
    hold on
    plot(simulatie.J_PLC_value(nummer,:),simulatie.opslagmaximum(nummer,:).^2 ./ (simulatie.opslagmaximum(nummer,:).^2 + 2^2), 'color',kleur(i))
    plot(simulatie.J_PLC_value(nummer,:),simulatie.opslagminimum(nummer,:).^2 ./ (simulatie.opslagminimum(nummer,:).^2 + 2^2), 'color',kleur(i))
    xlabel('J PLC [\muM/s]'); title('storagepart [-]')
    simulatie.description(nummer,:);
    if length(legendlist)==0;
        legendlist=[simulatie.description(nummer,:)];
    else
        legendlist=[legendlist ',' simulatie.description(nummer,:)];
    end
    
end
legenda = strsplit(legendlist,',');
legend(legenda);


%lijst=[2];
legendlist=[];
for i=[1:length(lijst)]
    nummer=lijst(i);
    figure(3)
    subplot(2,3,1)
    kleur=['k','b','r','g'];
    hold on
    plot(simulatie.J_PLC_value(nummer,:),simulatie.radiusmaximum(nummer,:)*1e6, 'color',kleur(i))
    plot(simulatie.J_PLC_value(nummer,:),simulatie.radiusminimum(nummer,:)*1e6, 'color',kleur(i))
    xlabel('J PLC [\muM/s]'); title('Radius [\mum]')
    simulatie.description(nummer,:);
    if length(legendlist)==0;
        legendlist=[simulatie.description(nummer,:)];
    else
        legendlist=[legendlist ',' simulatie.description(nummer,:)];
    end
    
end
legenda = strsplit(legendlist,',');
legend(legenda);

%lijst=[2];
legendlist=[];
for i=[1:length(lijst)]
    nummer=lijst(i);
    figure(3)
    subplot(2,3,2)
    kleur=['k','b','r','g'];
    hold on
    plot(simulatie.J_PLC_value(nummer,:),simulatie.NOmaximum(nummer,:), 'color',kleur(i))
    plot(simulatie.J_PLC_value(nummer,:),simulatie.NOminimum(nummer,:), 'color',kleur(i))
    xlabel('J PLC [\muM/s]'); title('NO concentration [\muM]')
    simulatie.description(nummer,:);
    if length(legendlist)==0;
        legendlist=[simulatie.description(nummer,:)];
    else
        legendlist=[legendlist ',' simulatie.description(nummer,:)];
    end
    
end
legenda = strsplit(legendlist,',');
legend(legenda);

%lijst=[2];
legendlist=[];
for i=[1:length(lijst)]
    nummer=lijst(i);
    figure(3)
    subplot(2,3,3)
    kleur=['k','b','r','g'];
    hold on
    plot(simulatie.J_PLC_value(nummer,:),simulatie.eNOSmaximum(nummer,:), 'color',kleur(i))
    plot(simulatie.J_PLC_value(nummer,:),simulatie.eNOSminimum(nummer,:), 'color',kleur(i))
    xlabel('J PLC [\muM/s]'); title('eNOS [\muM]')
    simulatie.description(nummer,:);
    if length(legendlist)==0;
        legendlist=[simulatie.description(nummer,:)];
    else
        legendlist=[legendlist ',' simulatie.description(nummer,:)];
    end
    
end
legenda = strsplit(legendlist,',');
legend(legenda);

%lijst=[2];
legendlist=[];
for i=[1:length(lijst)]
    nummer=lijst(i);
    figure(3)
    subplot(2,3,4)
    kleur=['k','b','r','g'];
    hold on
    plot(simulatie.J_PLC_value(nummer,:),simulatie.wssmaximum(nummer,:), 'color',kleur(i))
    plot(simulatie.J_PLC_value(nummer,:),simulatie.wssminimum(nummer,:), 'color',kleur(i))
    xlabel('J PLC [\muM/s]'); title('wss ')
    simulatie.description(nummer,:);
    if length(legendlist)==0;
        legendlist=[simulatie.description(nummer,:)];
    else
        legendlist=[legendlist ',' simulatie.description(nummer,:)];
    end
    
end
legenda = strsplit(legendlist,',');
legend(legenda);

%lijst=[2];
legendlist=[];
for i=[1:length(lijst)]
    nummer=lijst(i);
    figure(3)
    subplot(2,3,5)
    kleur=['k','b','r','g'];
    hold on
    plot(simulatie.J_PLC_value(nummer,:),simulatie.vimaximum(nummer,:), 'color',kleur(i))
    plot(simulatie.J_PLC_value(nummer,:),simulatie.viminimum(nummer,:), 'color',kleur(i))
    xlabel('J PLC [\muM/s]'); title('vi ')
    simulatie.description(nummer,:);
    if length(legendlist)==0;
        legendlist=[simulatie.description(nummer,:)];
    else
        legendlist=[legendlist ',' simulatie.description(nummer,:)];
    end
    
end
legenda = strsplit(legendlist,',');
legend(legenda);

%lijst=[2];
legendlist=[];
for i=[1:length(lijst)]
    nummer=lijst(i);
    figure(3)
    subplot(2,3,6)
    kleur=['k','b','r','g'];
    hold on
    plot(simulatie.J_PLC_value(nummer,:),simulatie.nNOSmaximum(nummer,:), 'color',kleur(i))
    plot(simulatie.J_PLC_value(nummer,:),simulatie.nNOSminimum(nummer,:), 'color',kleur(i))
    xlabel('J PLC [\muM/s]'); title('nNOS [\muM/s] ')
    simulatie.description(nummer,:);
    if length(legendlist)==0;
        legendlist=[simulatie.description(nummer,:)];
    else
        legendlist=[legendlist ',' simulatie.description(nummer,:)];
    end
    
end
legenda = strsplit(legendlist,',');
legend(legenda);

%lijst=[2];
legendlist=[];
for i=[1:length(lijst)]
    nummer=lijst(i);
    figure(4)
    subplot(2,3,1)
    kleur=['k','b','r','g'];
    hold on
    plot(simulatie.J_PLC_value(nummer,:),simulatie.JIP3maximum(nummer,:), 'color',kleur(i))
    plot(simulatie.J_PLC_value(nummer,:),simulatie.JIP3minimum(nummer,:), 'color',kleur(i))
    plot(simulatie.J_PLC_value(nummer,:),simulatie.JIP3mean(nummer,:), 'color',kleur(i))
    xlabel('J PLC [\muM/s]'); title('JIP3 [\muM/s]')
    simulatie.description(nummer,:);
    if length(legendlist)==0;
        legendlist=[simulatie.description(nummer,:)];
    else
        legendlist=[legendlist ',' simulatie.description(nummer,:)];
    end
    
end
legenda = strsplit(legendlist,',');
legend(legenda);

%lijst=[2];
legendlist=[];
for i=[1:length(lijst)]
    nummer=lijst(i);
    figure(4)
    subplot(2,3,2)
    kleur=['k','b','r','g'];
    hold on
    plot(simulatie.J_PLC_value(nummer,:),simulatie.Jextrusionmaximum(nummer,:), 'color',kleur(i))
    plot(simulatie.J_PLC_value(nummer,:),simulatie.Jextrusionminimum(nummer,:), 'color',kleur(i))
    plot(simulatie.J_PLC_value(nummer,:),simulatie.Jextrusionmean(nummer,:), 'color',kleur(i))
    xlabel('J PLC [\muM/s]'); title('J extrusion [\muM/s]')
    simulatie.description(nummer,:);
    if length(legendlist)==0;
        legendlist=[simulatie.description(nummer,:)];
    else
        legendlist=[legendlist ',' simulatie.description(nummer,:)];
    end
    
end
legenda = strsplit(legendlist,',');
legend(legenda);

%lijst=[2];
legendlist=[];
for i=[1:length(lijst)]
    nummer=lijst(i);
    figure(4)
    subplot(2,3,3)
    kleur=['k','b','r','g'];
    hold on
    plot(simulatie.J_PLC_value(nummer,:),simulatie.VOCCmaximum(nummer,:), 'color',kleur(i))
    plot(simulatie.J_PLC_value(nummer,:),simulatie.VOCCminimum(nummer,:), 'color',kleur(i))
    plot(simulatie.J_PLC_value(nummer,:),simulatie.VOCCmean(nummer,:), 'color',kleur(i))
    xlabel('J PLC [\muM/s]'); title('VOCC [\muM/s ')
    simulatie.description(nummer,:);
    if length(legendlist)==0;
        legendlist=[simulatie.description(nummer,:)];
    else
        legendlist=[legendlist ',' simulatie.description(nummer,:)];
    end
    
end
legenda = strsplit(legendlist,',');
legend(legenda);

%lijst=[2];
legendlist=[];
for i=[1:length(lijst)]
    nummer=lijst(i);
    figure(4)
    subplot(2,3,4)
    kleur=['k','b','r','g'];
    hold on
    plot(simulatie.J_PLC_value(nummer,:),simulatie.JNaCamaximum(nummer,:), 'color',kleur(i))
    plot(simulatie.J_PLC_value(nummer,:),simulatie.JNaCaminimum(nummer,:), 'color',kleur(i))
    plot(simulatie.J_PLC_value(nummer,:),simulatie.JNaCamean(nummer,:), 'color',kleur(i))
    xlabel('J PLC [\muM/s]'); title('JNaCa [\muM/s] ')
    simulatie.description(nummer,:);
    if length(legendlist)==0;
        legendlist=[simulatie.description(nummer,:)];
    else
        legendlist=[legendlist ',' simulatie.description(nummer,:)];
    end
    
end
legenda = strsplit(legendlist,',');
legend(legenda);



%lijst=[2];
legendlist=[];
for i=[1:length(lijst)]
    nummer=lijst(i);
    figure(4)
    subplot(2,3,5)
    kleur=['k','b','r','g'];
    hold on
    plot(simulatie.J_PLC_value(nummer,:),simulatie.stretchmaximum(nummer,:)*0.1, 'color',kleur(i))
    plot(simulatie.J_PLC_value(nummer,:),simulatie.stretchminimum(nummer,:)*0.1, 'color',kleur(i))
    plot(simulatie.J_PLC_value(nummer,:),simulatie.stretchmean(nummer,:)*0.1, 'color',kleur(i))
    xlabel('J PLC [\muM/s]'); title('stretch [\muM/s]')
    simulatie.description(nummer,:);
    if length(legendlist)==0;
        legendlist=[simulatie.description(nummer,:)];
    else
        legendlist=[legendlist ',' simulatie.description(nummer,:)];
    end
    
end
legenda = strsplit(legendlist,',');
legend(legenda);

%lijst=[2];
legendlist=[];
for i=[1:length(lijst)]
    nummer=lijst(i);
    figure(4)
    subplot(2,3,6)
    kleur=['k','b','r','g'];
    hold on
    plot(simulatie.J_PLC_value(nummer,:),simulatie.JCacouplingmaximum(nummer,:), 'color',kleur(i))
    plot(simulatie.J_PLC_value(nummer,:),simulatie.JCacouplingminimum(nummer,:), 'color',kleur(i))
    plot(simulatie.J_PLC_value(nummer,:),simulatie.JCacouplingmean(nummer,:), 'color',kleur(i))
    xlabel('J PLC [\muM/s]'); title('JCacoupling [\muM/s] ')
    simulatie.description(nummer,:);
    if length(legendlist)==0;
        legendlist=[simulatie.description(nummer,:)];
    else
        legendlist=[legendlist ',' simulatie.description(nummer,:)];
    end
    
end
legenda = strsplit(legendlist,',');
legend(legenda);

% 
% %lijst=[2];
% legendlist=[];
% for i=[1:length(lijst)]
%     nummer=lijst(i);
%     figure(5)
%     subplot(2,3,1)
%     kleur=['k','b','r','g'];
%     hold on
%     plot(simulatie.J_PLC_value(nummer,:),simulatie.radiusmaximum(nummer,:), 'color',kleur(i))
%     plot(simulatie.J_PLC_value(nummer,:),simulatie.radiusminimum(nummer,:), 'color',kleur(i))
%     plot(simulatie.J_PLC_value(nummer,:),simulatie.radiusmean(nummer,:), 'color',kleur(i))
%     xlabel('J PLC [\muM/s]'); title('radius [\mum]')
%     simulatie.description(nummer,:);
%     if length(legendlist)==0;
%         legendlist=[simulatie.description(nummer,:)];
%     else
%         legendlist=[legendlist ',' simulatie.description(nummer,:)];
%     end
%     
% end
% legenda = strsplit(legendlist,',');
% legend(legenda);
% 
% %lijst=[2];
% legendlist=[];
% for i=[1:length(lijst)]
%     nummer=lijst(i);
%     figure(5)
%     subplot(2,3,2)
%     kleur=['k','b','r','g'];
%     hold on
%     plot(simulatie.J_PLC_value(nummer,:),simulatie.Pressure(nummer,:)*10* 0.0075, 'color',kleur(i))
% %     plot(simulatie.J_PLC_value(nummer,:),simulatie.Jextrusionminimum(nummer,:), 'color',kleur(i))
% %     plot(simulatie.J_PLC_value(nummer,:),simulatie.Jextrusionmean(nummer,:), 'color',kleur(i))
%     xlabel('J PLC [\muM/s]'); title('sigma pr/h [Pa]')
%     simulatie.description(nummer,:);
%     if length(legendlist)==0;
%         legendlist=[simulatie.description(nummer,:)];
%     else
%         legendlist=[legendlist ',' simulatie.description(nummer,:)];
%     end
%     
% end
% legenda = strsplit(legendlist,',');
% legend(legenda);
% 
% %lijst=[2];
% legendlist=[];
% for i=[1:length(lijst)]
%     nummer=lijst(i);
%     figure(5)
%     subplot(2,3,3)
%     kleur=['k','b','r','g'];
%     hold on
%     plot(simulatie.J_PLC_value(nummer,:),exp(-1*7.4e-3*(simulatie.Pressure(nummer,:)*10* 0.0075-500)), 'color',kleur(i))
% %     plot(simulatie.J_PLC_value(nummer,:),simulatie.VOCCminimum(nummer,:), 'color',kleur(i))
% %     plot(simulatie.J_PLC_value(nummer,:),simulatie.VOCCmean(nummer,:), 'color',kleur(i))
%     xlabel('J PLC [\muM/s]'); title('e^(-alpha(sigma-sigmanul)) ')
%     simulatie.description(nummer,:);
%     if length(legendlist)==0;
%         legendlist=[simulatie.description(nummer,:)];
%     else
%         legendlist=[legendlist ',' simulatie.description(nummer,:)];
%     end
%     
% end
% legenda = strsplit(legendlist,',');
% legend(legenda);
% 
% %lijst=[2];
% legendlist=[];
% for i=[1:length(lijst)]
%     nummer=lijst(i);
%     figure(5)
%     subplot(2,3,4)
%     kleur=['k','b','r','g'];
%     hold on
%     plot(simulatie.J_PLC_value(nummer,:),6.1e-3./(1+ exp(-1*7.4e-3*(simulatie.Pressure(nummer,:)*10*0.0075-500))) , 'color',kleur(i))
% %     plot(simulatie.J_PLC_value(nummer,:),simulatie.JNaCaminimum(nummer,:), 'color',kleur(i))
% %     plot(simulatie.J_PLC_value(nummer,:),simulatie.JNaCamean(nummer,:), 'color',kleur(i))
%     xlabel('J PLC [\muM/s]'); title(' ')
%     simulatie.description(nummer,:);
%     if length(legendlist)==0;
%         legendlist=[simulatie.description(nummer,:)];
%     else
%         legendlist=[legendlist ',' simulatie.description(nummer,:)];
%     end
%     
% end
% legenda = strsplit(legendlist,',');
% legend(legenda);
% 
% %lijst=[2];
% legendlist=[];
% for i=[1:length(lijst)]
%     nummer=lijst(i);
%     figure(5)
%     subplot(2,3,5)
%     kleur=['k','b','r','g'];
%     hold on
%     plot(simulatie.J_PLC_value(nummer,:),simulatie.vimaximum(nummer,:), 'color',kleur(i))
%     plot(simulatie.J_PLC_value(nummer,:),simulatie.viminimum(nummer,:), 'color',kleur(i))
%     plot(simulatie.J_PLC_value(nummer,:),simulatie.vimean(nummer,:), 'color',kleur(i))
%     xlabel('J PLC [\muM/s]'); title('membrame potential [] ')
%     simulatie.description(nummer,:);
%     if length(legendlist)==0;
%         legendlist=[simulatie.description(nummer,:)];
%     else
%         legendlist=[legendlist ',' simulatie.description(nummer,:)];
%     end
%     
% end
% legenda = strsplit(legendlist,',');
% legend(legenda);
% 
% 
% %lijst=[2];
% legendlist=[];
% for i=[1:length(lijst)]
%     nummer=lijst(i);
%     figure(5)
%     subplot(2,3,6)
%     kleur=['k','b','r','g'];
%     hold on
%     plot(simulatie.J_PLC_value(nummer,:),simulatie.stretchmaximum(nummer,:)*0.1, 'color',kleur(i))
%     plot(simulatie.J_PLC_value(nummer,:),simulatie.stretchminimum(nummer,:)*0.1, 'color',kleur(i))
%     plot(simulatie.J_PLC_value(nummer,:),simulatie.stretchmean(nummer,:)*0.1, 'color',kleur(i))
%     xlabel('J PLC [\muM/s]'); title('stretch [\muM/s]')
%     simulatie.description(nummer,:);
%     if length(legendlist)==0;
%         legendlist=[simulatie.description(nummer,:)];
%     else
%         legendlist=[legendlist ',' simulatie.description(nummer,:)];
%     end
%     
% end
% legenda = strsplit(legendlist,',');
% legend(legenda);

