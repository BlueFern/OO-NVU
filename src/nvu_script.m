%% Demonstration script for new NVU model

%% Version control
% This is the latest version of the NVU model. It includes:
% - basic K+ / Ca2+ 

%% Construct NVU
% The NVU consists of a number of submodules, implemented as MATLAB
% classes, presently an astrocyte, a lumped SMC/EC model, and a model of
% the wall mechanics. 
%
% The parameters of each of these submodules are
% specified when the modules are constructed, here in the call to NVU, see
% for example the |SMCEC| part of the |NVU| call below:
%
% Options for the ODE solver (currently |ode15s|) are provided by
% specifying the |odeopts| parameter. The code works fine with default
% tolerances.
clear all
clc
odeopts = odeset('RelTol', 1e-03, 'AbsTol', 1e-03, 'MaxStep', 1, 'Vectorized', 1);

nv = NVU(Neuron('startpulse', 300, 'lengthpulse', 500), ...
    Astrocyte('startpulse', 300, 'lengthpulse', 500), ...
    WallMechanics(), ...
    SMCEC('J_PLC', 0.18), 'odeopts', odeopts);

nv.T = linspace(0, 1400, 1000);    
    
nv.simulate()

% Lists of quantities that can be retrieved from the NVU model after
% simulation are given in the documentation pages of the individual model
% components.
%
% Model components are retrieved using the |out| method as follows, where
% we do a simple plot of intracellular calcium from the SMC:

% Plot, e.g. Ca_i
% plot(nv.T, nv.out('Ca_i'))
% xlabel('time (s)')
% ylabel('[Ca^{2+}] (\muM)')

figure(2);
hold on
plot(nv.T, nv.out('R'))
xlabel('time (s)')
ylabel('R (m)')
hold off

figure(3);
hold on
plot(nv.T, nv.out('K_s'))
xlabel('time (s)')
ylabel('K_s')
hold off

%% Changing parameters, initial conditions, simulation time
% Changing parameters and inintial conditions involves adjusting properties
% of the model components of the |NVU| object.
%
% Parameters can be adjusted directly like so:
nv.smcec.params.J_PLC = 0.4;

%%
% Initial conditions can also be done in a single statement, but you end up
% with fairly long code lines, so the two-step method might be simpler if
% you're doing more than one:
nv.astrocyte.u0(nv.astrocyte.index.K_p) = 12;

% Alternative way -- has the same effect. These objects are passed by
% reference so you can do this
a = nv.astrocyte;
a.u0(a.index.K_p) = 12;

%%
% You can also adjust the simulation time:
nv.T = linspace(0, 1000, 2000);

%%
% To have make these changes have effect, you need to rerun the model
nv.simulate()

%% Recreate the wall plot from previous code
% We'll reset our parameters and time back to where they were. Don't forget
% you can look up default values and initial conditions on the model
% documentation pages.
nv.smcec.params.J_PLC = 0.18;
nv.astrocyte.u0(nv.astrocyte.index.K_p) = 3e3;
nv.T = linspace(0, 700, 1400);
nv.simulate();

subplot(3, 2, 1)
plot(nv.T, nv.out('M'))
xlabel('Time')
ylabel('Fraction [-]')
title('[M]')

subplot(3, 2, 2)
plot(nv.T, nv.out('Mp'))
xlabel('Time')
ylabel('Fraction [-]')
title('[Mp]')

subplot(3, 2, 3)
plot(nv.T, nv.out('AMp'))
xlabel('Time')
ylabel('Fraction [-]')
title('[AMp]')

subplot(3, 2, 4)
plot(nv.T, nv.out('AM'))
xlabel('Time')
ylabel('Fraction [-]')
title('[AM]')

subplot(3, 2, 5)
plot(nv.T, nv.out('F_r'))
xlabel('Time')
ylabel('Fraction [-]')
title('F_r')

subplot(3, 2, 6)
plot(nv.T, 1e6 * nv.out('R'))
xlabel('Time')
ylabel('\mu m')
title('Radius')

%% Disable one of the variables
% It is straightforward to disable one of the variables. This sets its ODE
% to du/dt = 0, so it remains at its initial condition for the duration of
% the simulation.
nv.wall.enabled(nv.wall.index.R) = false;
nv.simulate()
subplot(3, 2, 1)
plot(nv.T, nv.out('M'))
xlabel('Time')
ylabel('Fraction [-]')
title('[M]')

subplot(3, 2, 2)
plot(nv.T, nv.out('Mp'))
xlabel('Time')
ylabel('Fraction [-]')
title('[Mp]')

subplot(3, 2, 3)
plot(nv.T, nv.out('AMp'))
xlabel('Time')
ylabel('Fraction [-]')
title('[AMp]')

subplot(3, 2, 4)
plot(nv.T, nv.out('AM'))
xlabel('Time')
ylabel('Fraction [-]')
title('[AM]')

subplot(3, 2, 5)
plot(nv.T, nv.out('F_r'))
xlabel('Time')
ylabel('Fraction [-]')
title('F_r')

subplot(3, 2, 6)
plot(nv.T, 1e6 * nv.out('R'))
xlabel('Time')
ylabel('\mu m')
title('Radius')

%% Switch TRPV4 channel off:
nv.astrocyte.params.switchBK = 1;
nv.astrocyte.params.reverseBK = 0; % V
nv.astrocyte.params.G_BK_k = 4.3e3; % pS (later converted to mho m^-2)
nv.astrocyte.params.trpv_switch = 0;
nv.astrocyte.params.epshalf_k = 0.16;
nv.astrocyte.params.Ca_4 = 0.15; % uM
nv.astrocyte.params.v_7 = -15e-3; % V

% Re-run simulation:
nv.simulate()

%% Switch TRPV4 channel back on (default): 
% This set of parameters is chosen to be the "correct" one. For variations,
% see nvu_script_easy_run by Joerik!
nv.astrocyte.params.switchBK = 0;
nv.astrocyte.params.reverseBK = -0.08135; % V
nv.astrocyte.params.G_BK_k = 225; % pS (later converted to mho m^-2)
nv.astrocyte.params.trpv_switch = 1;
nv.astrocyte.params.epshalf_k = 0.1;
nv.astrocyte.params.Ca_4 = 0.35; % uM
nv.astrocyte.params.v_7 = -13.57e-3; % V

% Re-run simulation:
nv.simulate()

%% plot all state variables:
for i = 1:2
    figure(1)
    set(gcf,'name','Neuron and Astrocyte')
    subplot(5,4,i)
    plot(nv.T, nv.out(nv.neuron.varnames{i}))
    h1(i) = gca();
    ylabel(nv.neuron.varnames{i})
    hold on
end
    linkaxes(h1, 'x');

for i = 1:17
    figure(1)
    subplot(5,4,2+i)
    plot(nv.T, nv.out(nv.astrocyte.varnames{i}))
    h1(i) = gca();
    ylabel(nv.astrocyte.varnames{i})
    hold on
end
    linkaxes(h1, 'x');


for i = 1:10
    figure(2)
    set(gcf,'name','ECSMC & WallMechanics')
    subplot(3,5,i)
    plot(nv.T, nv.out(nv.smcec.varnames{i}))
    h2(i) = gca();
    ylabel(nv.smcec.varnames{i})
    hold on
end

    
for i = 1:4
    figure(2)
    subplot(3,5,10+i)
    plot(nv.T, nv.out(nv.wall.varnames{i}))
    h2(10+i) = gca();
    ylabel(nv.wall.varnames{i})
    hold on
end
    linkaxes(h2, 'x');
    hold on
