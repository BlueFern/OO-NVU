%% Loop made for creating a figure that shows calcium concentration over 
%% increasing J_PLC for new NVU model version 1.2

simulatie.nummer=[6 7 8 9];%
simulatie.description   = ['old - normal              ';     %1
                           'new - normal              ';     %2
                           'new - VOCC*0.5            ';     %3
                           'new - VOCC*0.60           ';     %4
                           'new - VOCC*0.75           ';     %5
                           'new - VOCC*0.90           ';     %6
                           'new - VOCC*0.10           ';     %7
                           'new - VOCC*0.01           ';     %8
                           'new - PMCA*0.50           ';     %9
                           'new - PMCA*0.10           ';     %10
                           'new - PMCA*0.05           ';     %11
                           'new - PMCA*0.001          ';     %12
                           'new - PMCA*0.005          ';     %13
                           'old - extr*0.50           ';     %14
                           'old - extr*0.10           ';     %15
                           'old - extr*0.05           ';     %16
                           'old - extr*0.001          ';     %17
                           'old - extr*0.005          ';     %18
                           'old - VOCC*0.5            ';     %19
                           'old - VOCC*0.60           ';     %20
                           'old - VOCC*0.75           ';     %21
                           'old - VOCC*0.90           ';     %22
                           'old - VOCC*0.10           ';     %23
                           'old - VOCC*0.01           ';     %24
                           'old - no neural           ';     %25
                           'new - no neural           ';     %26
                           'new - VOCC*0.10  PMCA*0.10';     %27
                           'new - pressure*0.75       ';     %28
                           'new - pressure*0.50       ';     %29
                           'new - pressure*0.25       ';     %30
                           'new - pressure*0.10       ';     %31
                           'new - pressure*0.05       ';     %32
                           'new - CICR*0.5            ';     %33
                           'new - CICR*0.1            ';     %34
                           'new - CICR*0.05           ';     %35                                    
                           'old - pressure*0.75       ';     %36
                           'old - pressure*0.50       ';     %37
                           'old - pressure*0.25       ';     %38
                           'old - pressure*0.10       ';     %39
                           'old - pressure*2          ';     %40
                           'old - extr*0.50- no neural';     %41
                           'old - extr*0.10- no neural';     %42
                           'old - extr*0.05- no neural';     %43
                           'old - extr*0.001-no neural';     %44
                           'old - extr*0.005-no neural'];    %45                       
simulatie.settings      = ['nv.SMCEC.params.newextrusionswitch';
                           'nv.SMCEC.params.G_Ca_i            ';
                           'nv.SMCEC.params.I_PMCA            ';
                           'nv.SMCEC.params.D_i               ';
                           'startpulse                        ';
                           'pressure                          ';
                           'CICR                              '];
 VOCC=1.29e-3;
 IPMCA=5.37;
 D_i=0.24;
 normalstart=100;
 simpres=4000;
 CICR=55;
simulatie.settingsvalue.newextrusion = [0 1 1 1 1 1 1 1 1 1 1 1 1 ...
                                        0 0 0 0 0 0 0 0 0 0 0 0 1 ...
                                        1 1 1 1 1 1 1 1 1 0 0 0 0 ...
                                        0 1 1 1 1 1];
simulatie.settingsvalue.G_Ca_i       = VOCC*[1 1 0.5 0.6 0.75 0.90 0.10 0.01 1 1 1 1 1 ...
                                             1 1 1 1 1 0.5 0.6 0.75 0.90 0.10 0.01 1 1 ...
                                             0.10 1 1 1 1 1 1 1 1 1 1 1 1 ...
                                             1 1 1 1 1 1];
simulatie.settingsvalue.I_PMCA       = IPMCA*[1 1 1 1 1 1 1 1 0.5 0.1 0.05 0.001 0.005 ...
                                              1 1 1 1 1 1 1 1 1 1 1 1 1 ...
                                              0.10 1 1 1 1 1 1 1 1 1 1 1 1 ...
                                              1 1 1 1 1 1];
simulatie.settingsvalue.D_i          = D_i*[1 1 1 1 1 1 1 1 1 1 1 1 1 ...
                                            0.50 0.10 0.05 0.001 0.005 1 1 1 1 1 1 1 1 ...
                                            1 1 1 1 1 1 1 1 1 1 1 1 1 ...
                                            1 0.50 0.10 0.05 0.001 0.005];
simulatie.settingsvalue.startpulse   = normalstart+[0 0 0 0 0 0 0 0 0 0 0 0 0 ...
                                                    0 0 0 0 0 0 0 0 0 0 0 900 900 ...
                                                    0 0 0 0 0 0 0 0 0 1 1 1 1 ...
                                                    1 900 900 900 900 900];
simulatie.settingsvalue.pressure     = simpres*[1 1 1 1 1 1 1 1 1 1 1 1 1 ...
                                                1 1 1 1 1 1 1 1 1 1 1 1 1 ...
                                                1 0.75 0.50 0.25 0.10 0.05 1 1 1 0.75 0.50 0.25 0.10 ...
                                                2 1 1 1 1 1];
simulatie.settingsvalue.CICR         =    CICR*[1 1 1 1 1 1 1 1 1 1 1 1 1 ...
                                                1 1 1 1 1 1 1 1 1 1 1 1 1 ...
                                                1 1 1 1 1 1 0.5 0.1 0.05 1 1 1 1 ...
                                                1 1 1 1 1 1];
for xxx=simulatie.nummer
    eenmaximumlist=[];
    eenminimumlist=[];
    meanlist=[];
    opslagmaximum=[];
    opslagminimum=[];
    opslagmean=[];
    radiusmaximum=[];
    radiusminimum=[];
    radiusmean=[];
    NOmaximum=[];
    NOminimum=[];
    NOmean=[];
    eNOSmaximum=[];
    eNOSminimum=[];
    eNOSmean=[];
    nNOSmaximum=[];
    nNOSminimum=[];
    nNOSmean=[];    
    wssmaximum=[];
    wssminimum=[];
    wssmean=[];
    vimaximum=[];
    viminimum=[];
    vimean=[];
    VOCCmaximum=[];
    VOCCminimum=[];
    VOCCmean=[];
    stretchmaximum=[];
    stretchminimum=[];
    stretchmean=[];
    Jextrusionmaximum=[];
    Jextrusionminimum=[];
    Jextrusionmean=[];
    JNaCamaximum=[];
    JNaCaminimum=[];
    JNaCamean=[];
    JCacouplingmaximum=[];
    JCacouplingminimum=[];
    JCacouplingmean=[];
    JIP3maximum=[];
    JIP3minimum=[];
    JIP3mean=[];
    
    J_PLC_value=[0.1:0.05:0.8];
    %simulatie.settings(1,:)=simulatie.settingsvalue.newextrusion(xxx);
    %simulatie.settings(2,:)=simulatie.settingsvalue.G_Ca_i(xxx);
    Calciumconcentration=[];
    for xx=J_PLC_value
    
        % Things to change:
        %XLIM1 = 50; XLIM2 = 500;            % x limits of any plots
        PRESSURE = simulatie.settingsvalue.pressure(xxx);% Pressure in Pa, change here!! Because the parameter is used in both SMCEC and WallMechanics ***********, 4000 Pa is default
        PRESSURE_CHANGE=2000;                % Pressure at 100 sec till 300 seconds when switch is on, 2000 Pa is default
        START_PULSE  = simulatie.settingsvalue.startpulse(xxx);                  % Start time of neuronal stimulation
        LENGTH_PULSE = 1900;       	         % Length of stimulation 
        GLU_SWITCH   = 0;                    % Turn on glutamate input to SC, default 1
        NO_SWITCH    = 1;                    % Turn on Nitric Oxide production, default 1
        TRPV4_SWITCH = 1;                    % Turn on TRPV4 channel flux in AC, default 1
        PLC_SWITCH   = 0;                    % Turn on PLC change in time, default 1
        PressureSwitch =0;                   % Turn on pressure change in time, default 1
        J_PLC_steadystate = xx;    % Jplc value in EC: 0.18 for steady state (default)
        J_PLC_vasomotion  = xx;     % Jplc value in EC:  0.4 for vasomotion   (default) from 100 to 300 seconds
        NEWEXTRUSIONSWITCH = simulatie.settingsvalue.newextrusion(xxx);      %Turn on new extrusion based on Kapela, default 1
        
        odeopts = odeset('RelTol', 1e-04, 'AbsTol', 1e-04, 'MaxStep', 0.5, 'Vectorized', 1);


        % Initiate the NVU
        nv = NVU(Neuron('startpulse', START_PULSE, 'lengthpulse', LENGTH_PULSE, 'GluSwitch', GLU_SWITCH, 'NOswitch', NO_SWITCH), ...
            Astrocyte('startpulse', START_PULSE, 'lengthpulse', LENGTH_PULSE, 'rhoSwitch', GLU_SWITCH, 'TRPV4switch', TRPV4_SWITCH), ...
            WallMechanics('PressureSwitch', PressureSwitch,'P_T', PRESSURE,'Pressure_change', PRESSURE_CHANGE), ...
            SMCEC('PressureSwitch', PressureSwitch, 'J_PLC_steadystate', J_PLC_steadystate, 'J_PLC_vasomotion',J_PLC_vasomotion, 'newextrusionswitch', NEWEXTRUSIONSWITCH,  ...
            'G_Ca_i',simulatie.settingsvalue.G_Ca_i(xxx),  'I_PMCA', simulatie.settingsvalue.I_PMCA(xxx), 'NOswitch', NO_SWITCH, 'P_T', PRESSURE, 'PLCSwitch', PLC_SWITCH, ...
            'Pressure_change', PRESSURE_CHANGE, 'D_i',simulatie.settingsvalue.D_i(xxx), 'C_i',simulatie.settingsvalue.CICR(xxx) ), 'odeopts', odeopts);


        nv.T = linspace(0, 1000, 10000);     % Initiate time vector
        nv.simulate(); % Run simulation


        Calciumconcentration=nv.out('Ca_i');
        eenmaximumlist=[eenmaximumlist max(Calciumconcentration(9000:9950))];
        eenminimumlist=[eenminimumlist min(Calciumconcentration(9000:9950))];
        meanlist=[meanlist mean(Calciumconcentration(9000:9950))];
        opslag=nv.out('s_i');
        opslagmaximum=[opslagmaximum max(opslag(9000:9950))];
        opslagminimum=[opslagminimum min(opslag(9000:9950))];
        opslagmean=[opslagmean mean(opslag(9000:9950))];
        radius=nv.out('R');
        radiusmaximum=[radiusmaximum max(radius(9000:9950))];
        radiusminimum=[radiusminimum min(radius(9000:9950))];
        radiusmean=[radiusmean mean(radius(9000:9950))];
        NO=nv.out('NO_i');
        NOmaximum=[NOmaximum max(NO(9000:9950))];
        NOminimum=[NOminimum min(NO(9000:9950))];
        NOmean=[NOmean mean(NO(9000:9950))];
        eNOS=nv.out('eNOS_act_j');
        eNOSmaximum=[eNOSmaximum max(eNOS(9000:9950))];
        eNOSminimum=[eNOSminimum min(eNOS(9000:9950))];
        eNOSmean=[eNOSmean mean(eNOS(9000:9950))];
        nNOS=nv.out('eNOS_act_j');
        nNOSmaximum=[nNOSmaximum max(nNOS(9000:9950))];
        nNOSminimum=[nNOSminimum min(nNOS(9000:9950))];
        nNOSmean=[nNOSmean mean(nNOS(9000:9950))];        
        wss=nv.out('tau_wss');
        wssmaximum=[wssmaximum max(wss(9000:9950))];
        wssminimum=[wssminimum min(wss(9000:9950))];
        wssmean=[wssmean mean(wss(9000:9950))];
        vi=nv.out('v_i');
        vimaximum=[vimaximum max(vi(9000:9950))];
        viminimum=[viminimum min(vi(9000:9950))];
        vimean=[vimean mean(vi(9000:9950))];
        VOCC=nv.out('J_VOCC_i');
        VOCCmaximum=[VOCCmaximum max(VOCC(9000:9950))];
        VOCCminimum=[VOCCminimum min(VOCC(9000:9950))];
        VOCCmean=[VOCCmean mean(VOCC(9000:9950))];
        stretch=nv.out('J_stretch_i');
        stretchmaximum=[stretchmaximum max(stretch(9000:9950))];
        stretchminimum=[stretchminimum min(stretch(9000:9950))];
        stretchmean=[stretchmean mean(stretch(9000:9950))];

        extrusion=nv.out('J_extrusion_i');
        Jextrusionmaximum=[Jextrusionmaximum max(extrusion(9000:9950))];
        Jextrusionminimum=[Jextrusionminimum min(extrusion(9000:9950))];
        Jextrusionmean=[Jextrusionmean mean(extrusion(9000:9950))];
        
        NaCa=nv.out('J_NaCa_i');
        JNaCamaximum=[JNaCamaximum max(NaCa(9000:9950))];
        JNaCaminimum=[JNaCaminimum min(NaCa(9000:9950))];
        JNaCamean=[JNaCamean mean(NaCa(9000:9950))];        
        
        Cacoupling=nv.out('J_Ca_coup_i');
        JCacouplingmaximum=[JCacouplingmaximum max(Cacoupling(9000:9950))];
        JCacouplingminimum=[JCacouplingminimum min(Cacoupling(9000:9950))];
        JCacouplingmean=[JCacouplingmean mean(Cacoupling(9000:9950))];
        
        IP3=nv.out('J_IP3_i');
        JIP3maximum=[JIP3maximum max(IP3(9000:9950))];
        JIP3minimum=[JIP3minimum min(IP3(9000:9950))];
        JIP3mean=[JIP3mean mean(IP3(9000:9950))];        
        
    end
    simulatie.maximumlist(xxx,:)=eenmaximumlist;
    simulatie.minimumlist(xxx,:)=eenminimumlist;
    simulatie.meanlist(xxx,:)=meanlist;
    simulatie.J_PLC_value(xxx,:)=J_PLC_value;
    simulatie.opslagmaximum(xxx,:)=opslagmaximum;
    simulatie.opslagminimum(xxx,:)=opslagminimum;
    simulatie.opslagmean(xxx,:)=opslagmean;
    simulatie.radiusmaximum(xxx,:)=radiusmaximum;
    simulatie.radiusminimum(xxx,:)=radiusminimum;
    simulatie.radiusmean(xxx,:)=radiusmean;
    simulatie.NOmaximum(xxx,:)=NOmaximum;
    simulatie.NOminimum(xxx,:)=NOminimum;
    simulatie.NOmean(xxx,:)=NOmean;
    simulatie.eNOSmaximum(xxx,:)=eNOSmaximum;
    simulatie.eNOSminimum(xxx,:)=eNOSminimum;
    simulatie.eNOSmean(xxx,:)=eNOSmean;
    simulatie.nNOSmaximum(xxx,:)=nNOSmaximum;
    simulatie.nNOSminimum(xxx,:)=nNOSminimum;
    simulatie.nNOSmean(xxx,:)=nNOSmean;
    simulatie.wssmaximum(xxx,:)=wssmaximum;
    simulatie.wssminimum(xxx,:)=wssminimum;
    simulatie.wssmean(xxx,:)=wssmean;
    simulatie.vimaximum(xxx,:)=vimaximum;
    simulatie.viminimum(xxx,:)=viminimum;
    simulatie.vimean(xxx,:)=vimean;    
    simulatie.VOCCmaximum(xxx,:)=VOCCmaximum;
    simulatie.VOCCminimum(xxx,:)=VOCCminimum;
    simulatie.VOCCmean(xxx,:)=VOCCmean;
    simulatie.stretchmaximum(xxx,:)=stretchmaximum;
    simulatie.stretchminimum(xxx,:)=stretchminimum;
    simulatie.stretchmean(xxx,:)=stretchmean;
    
    simulatie.Jextrusionmaximum(xxx,:)=Jextrusionmaximum;
    simulatie.Jextrusionminimum(xxx,:)=Jextrusionminimum;
    simulatie.Jextrusionmean(xxx,:)=Jextrusionmean;   
    
    simulatie.JNaCamaximum(xxx,:)=JNaCamaximum;
    simulatie.JNaCaminimum(xxx,:)=JNaCaminimum;
    simulatie.JNaCamean(xxx,:)=JNaCamean;
    
    simulatie.JCacouplingmaximum(xxx,:)=JCacouplingmaximum;
    simulatie.JCacouplingminimum(xxx,:)=JCacouplingminimum;
    simulatie.JCacouplingmean(xxx,:)=JCacouplingmean;
    
    simulatie.JIP3maximum(xxx,:)=JIP3maximum;
    simulatie.JIP3minimum(xxx,:)=JIP3minimum;
    simulatie.JIP3mean(xxx,:)=JIP3mean;
    
end

% save('simulatiemeer_metPa_2','simulatie')


% 
% 
% lijst=[1 36 37];
% legendlist=[];
% for i=[1:length(lijst)]
%     nummer=lijst(i);
%     figure(1)
%     kleur=['k','b','r','g'];
%     hold on
%     plot(simulatie.J_PLC_value(nummer,:),simulatie.maximumlist(nummer,:), 'color',kleur(i))
%     plot(simulatie.J_PLC_value(nummer,:),simulatie.minimumlist(nummer,:), 'color',kleur(i))
%     xlabel('J PLC [\muM/s]'); title('Ca2+ concentration SMC [\muM]')
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
%     figure(2)
%     subplot(2,1,1)
%     kleur=['k','b','r','g'];
%     hold on
%     plot(simulatie.J_PLC_value(nummer,:),simulatie.maximumlist(nummer,:).^4 ./ (simulatie.maximumlist(nummer,:).^4 + 0.9^4), 'color',kleur(i))
%     plot(simulatie.J_PLC_value(nummer,:),simulatie.minimumlist(nummer,:).^4 ./ (simulatie.minimumlist(nummer,:).^4 + 0.9^4), 'color',kleur(i))
%     xlabel('J PLC [\muM/s]'); title('calciumpart [-]')
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
%     figure(2)
%     subplot(2,1,2)
%     kleur=['k','b','r','g'];
%     hold on
%     plot(simulatie.J_PLC_value(nummer,:),simulatie.opslagmaximum(nummer,:).^2 ./ (simulatie.opslagmaximum(nummer,:).^2 + 2^2), 'color',kleur(i))
%     plot(simulatie.J_PLC_value(nummer,:),simulatie.opslagminimum(nummer,:).^2 ./ (simulatie.opslagminimum(nummer,:).^2 + 2^2), 'color',kleur(i))
%     xlabel('J PLC [\muM/s]'); title('storagepart [-]')
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
%     figure(3)
%     subplot(2,3,1)
%     kleur=['k','b','r','g'];
%     hold on
%     plot(simulatie.J_PLC_value(nummer,:),simulatie.radiusmaximum(nummer,:), 'color',kleur(i))
%     plot(simulatie.J_PLC_value(nummer,:),simulatie.radiusminimum(nummer,:), 'color',kleur(i))
%     xlabel('J PLC [\muM/s]'); title('Radius [\mum]')
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
%     figure(3)
%     subplot(2,3,2)
%     kleur=['k','b','r','g'];
%     hold on
%     plot(simulatie.J_PLC_value(nummer,:),simulatie.NOmaximum(nummer,:), 'color',kleur(i))
%     plot(simulatie.J_PLC_value(nummer,:),simulatie.NOminimum(nummer,:), 'color',kleur(i))
%     xlabel('J PLC [\muM/s]'); title('NO concentration [\muM]')
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
%     figure(3)
%     subplot(2,3,3)
%     kleur=['k','b','r','g'];
%     hold on
%     plot(simulatie.J_PLC_value(nummer,:),simulatie.eNOSmaximum(nummer,:), 'color',kleur(i))
%     plot(simulatie.J_PLC_value(nummer,:),simulatie.eNOSminimum(nummer,:), 'color',kleur(i))
%     xlabel('J PLC [\muM/s]'); title('eNOS concentration [\muM]')
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
%     figure(3)
%     subplot(2,3,4)
%     kleur=['k','b','r','g'];
%     hold on
%     plot(simulatie.J_PLC_value(nummer,:),simulatie.wssmaximum(nummer,:), 'color',kleur(i))
%     plot(simulatie.J_PLC_value(nummer,:),simulatie.wssminimum(nummer,:), 'color',kleur(i))
%     xlabel('J PLC [\muM/s]'); title('wss ')
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
%     figure(3)
%     subplot(2,3,5)
%     kleur=['k','b','r','g'];
%     hold on
%     plot(simulatie.J_PLC_value(nummer,:),simulatie.vimaximum(nummer,:), 'color',kleur(i))
%     plot(simulatie.J_PLC_value(nummer,:),simulatie.viminimum(nummer,:), 'color',kleur(i))
%     xlabel('J PLC [\muM/s]'); title('vi ')
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
%     figure(3)
%     subplot(2,3,6)
%     kleur=['k','b','r','g'];
%     hold on
%     plot(simulatie.J_PLC_value(nummer,:),simulatie.VOCCmaximum(nummer,:), 'color',kleur(i))
%     plot(simulatie.J_PLC_value(nummer,:),simulatie.VOCCminimum(nummer,:), 'color',kleur(i))
%     xlabel('J PLC [\muM/s]'); title('VOCC ')
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
