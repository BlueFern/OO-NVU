v_sa = linspace(-100, 20, 1000);

v_switch = -65;
slope = 0.1;
Glu_max = 1846;

Glu = 0.5*Glu_max*(1 + tanh((v_sa-v_switch)/slope));

figure(99);
plot(v_sa, Glu);