%% Loop made for creating a figure that shows calcium concentration over 
%% increasing J_PLC for new NVU model version 1.2
% file used to make a graph of changing parameters with increasing Jplc 
% 
simulatie.nummer=[1 3 5];%
simulatie.description   = ['old - normal                                ';     %1
                           'new - normal                                ';     %2                                   

                           'old - pressure*0.50                         ';     %3
                           'old - pressure*0.25                         ';     %4
                           'old - pressure*2                            ';     %5
                           
                           'old - normal nocw                           ';     %6
                           'old - pressure*0.50 nocw                    ';     %7
                           'old - pressure*0.25 nocw                    ';     %8
                           'old - pressure*2    nocw                    ';     %9
                           
                           'old - normal no NO                          ';     %10
                           'old - pressure*0.50 no NO                   ';     %11
                           'old - pressure*0.25 no NO                   ';     %12
                           'old - pressure*2    no NO                   ';     %13
                           
                           'old - normal no stretch                     ';     %14
                           'old - pressure*0.50 no stretch              ';     %15
                           'old - pressure*0.25 no stretch              ';     %16
                           'old - pressure*2    no stretch              ';     %17                          

                           'old - normal no TRPV                        ';     %18
                           'old - pressure*0.50 no TRPV                 ';     %19
                           'old - pressure*0.25 no TRPV                 ';     %20
                           'old - pressure*2    no TRPV                 ';     %21
                           
                           'old - normal no NO no stretch               ';     %22
                           'old - pressure*0.50 no NO no stretch        ';     %23
                           'old - pressure*0.25 no NO no stretch        ';     %24
                           'old - pressure*2 no NO no stretch           ';     %25
 
                           'old - normal no TRPV no stretch             ';     %26
                           'old - pressure*0.50 no TRPV no stretch      ';     %27
                           'old - pressure*0.25 no TRPV no stretch      ';     %28
                           'old - pressure*2 no TRPV no stretch         ';     %29
                           
                           'old - normal no TRPV no NO                  ';     %30
                           'old - pressure*0.50 no NO                   ';     %31
                           'old - pressure*0.25 no NO                   ';     %32   
                           'old - pressure*2 no TRPV NO                 ';     %33
                           
                           'old - normal no TRPV no NO no stretch       ';     %34
                           'old - pressure*0.50 no NO no NO no stretch  ';     %35
                           'old - pressure*0.25 no NO no NO no stretch  ';     %36
                           'old - pressure*2 no TRPV NO no NO no stretch';     %37
                           
                           'old - normal sigma*1.1                      ';     %38
                           'old - normal sigma*1.25                     ';     %39
                           'old - normal sigma*1.5                      ';     %40
                           'old - normal sigma*0.9                      ';     %41
                           'old - normal sigma*0.75                     ';     %42
                           'old - normal sigma*0.5                      ';     %43
                           
                           'old - pressure*0.50 sigma*1.1               ';     %44
                           'old - pressure*0.50 sigma*1.25              ';     %45
                           'old - pressure*0.50 sigma*1.5               ';     %46
                           'old - pressure*0.50 sigma*0.9               ';     %47
                           'old - pressure*0.50 sigma*0.75              ';     %48
                           'old - pressure*0.50 sigma*0.5               ';     %49                           
                           
                           'old - pressure*2 sigma*1.1                  ';     %50
                           'old - pressure*2 sigma*1.25                 ';     %51
                           'old - pressure*2 sigma*1.5                  ';     %52
                           'old - pressure*2 sigma*0.9                  ';     %53
                           'old - pressure*2 sigma*0.75                 ';     %54
                           'old - pressure*2 sigma*0.5                  '];    %55                         
                           
                           
                           
                           

 normalstart=100;
 simpres=4000;
 

simulatie.settingsvalue.startpulse   = normalstart+[0 0 0 0  ...
                                                    0 0 0 0  ...
                                                    0 0 0 0  ...
                                                    0 0 0 0  ...
                                                    0 0 0 0  ...
                                                    0 0 0 0  ...
                                                    0 0 0 0  ...
                                                    0 0 0 0  ...
                                                    0 0 0 0  ...
                                                    0 0 0 0  ...
                                                    0 0 0 0  ...
                                                    0 0 0 0  ...
                                                    0 0 0 0  ...
                                                    0 0 0];
                                                
simulatie.settingsvalue.pressure     = simpres*[1 1 0.50 0.25  ...
                                                2 1 0.50 0.25  ...
                                                2 1 0.50 0.25  ...
                                                2 1 0.50 0.25  ...
                                                2 1 0.50 0.25  ...
                                                2 1 0.50 0.25  ...
                                                2 1 0.50 0.25  ...
                                                2 1 0.50 0.25  ...
                                                2 1 0.50 0.25  ...
                                                2 1 1    1     ...
                                                1 1 1    0.50  ...
                                                0.50 0.50 0.50 0.50 ...
                                                0.50 2    2    2    ...
                                                2 2  2];


simulatie.settingsvalue.stretch      = [1 1 1 1  ...
                                        1 1 1 1  ...
                                        1 1 1 1  ...
                                        1 0 0 0  ...
                                        0 1 1 1  ...
                                        1 0 0 0  ...
                                        0 0 0 0  ...
                                        0 1 1 1  ...
                                        1 0 0 0  ...
                                        0 1 1 1  ...
                                        1 1 1 1  ...
                                        1 1 1 1  ...
                                        1 1 1 1  ...
                                        1 1 1];                           
                                            
                                              
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
    JCICRimaximum=[];
    JCICRiminimum=[];
    JCICRimean=[]; 
    Frmaximum=[];
    Frminimum=[];
    Frmean=[];
    nummerbezig=simulatie.description(xxx,:)
    J_PLC_value=[0.1:0.0025:0.8];
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
        GLU_SWITCH   = 1;                    % Turn on glutamate input to SC, default 1
        NO_SWITCH    = 1;                     % Turn on Nitric Oxide production, default 1
        TRPV4_SWITCH = 1;                     % Turn on TRPV4 channel flux in AC, default 1
        PLC_SWITCH   = 0;                    % Turn on PLC change in time, default 1
        PressureSwitch =0;                   % Turn on pressure change in time, default 1
        J_PLC_steadystate = xx;    % Jplc value in EC: 0.18 for steady state (default)
        J_PLC_vasomotion  = xx;     % Jplc value in EC:  0.4 for vasomotion   (default) from 100 to 300 seconds

        odeopts = odeset('RelTol', 1e-04, 'AbsTol', 1e-04, 'MaxStep', 0.5, 'Vectorized', 1);


        % Initiate the NVU
        nv = NVU(Neuron('startpulse', START_PULSE, 'lengthpulse', LENGTH_PULSE, 'GluSwitch', GLU_SWITCH, 'NOswitch', NO_SWITCH), ...
            Astrocyte('startpulse', START_PULSE, 'lengthpulse', LENGTH_PULSE, 'rhoSwitch', GLU_SWITCH, 'TRPV4switch', TRPV4_SWITCH), ...
            WallMechanics('PressureSwitch', PressureSwitch,'P_T', PRESSURE,'Pressure_change', PRESSURE_CHANGE), ...
            SMCEC('PressureSwitch', PressureSwitch, 'J_PLC_steadystate', J_PLC_steadystate, 'J_PLC_vasomotion',J_PLC_vasomotion,  ...
             'NOswitch', NO_SWITCH, 'P_T', PRESSURE, 'PLCSwitch', PLC_SWITCH, ...
            'Pressure_change', PRESSURE_CHANGE), 'odeopts', odeopts);


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
        nNOS=nv.out('nNOS_act_n');
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
        
        CICRi=nv.out('J_CICR_i');
        JCICRimaximum=[JCICRimaximum max(CICRi(9000:9950))];
        JCICRiminimum=[JCICRiminimum min(CICRi(9000:9950))];
        JCICRimean=[JCICRimean mean(CICRi(9000:9950))]; 
        
        Fr=nv.out('F_r');
        Frmaximum=[Frmaximum max(Fr(9000:9950))];
        Frminimum=[Frminimum min(Fr(9000:9950))];
        Frmean=[Frmean mean(Fr(9000:9950))];
        
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
    
    simulatie.Frmaximum(xxx,:)=Frmaximum;
    simulatie.Frminimum(xxx,:)=Frminimum;
    simulatie.Frmean(xxx,:)=Frmean;
    
    simulatie.JCICRimaximum(xxx,:)=JCICRimaximum;
    simulatie.JCICRiminimum(xxx,:)=JCICRiminimum;
    simulatie.JCICRimean(xxx,:)=JCICRimean;
    simulatie.Pressure(xxx,:)=simulatie.settingsvalue.pressure(xxx)*ones(1,length(J_PLC_value));
end
%save('simulatie_metPa_old_stretchandersom_10jan','simulatie')
