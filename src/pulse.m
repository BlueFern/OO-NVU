clear all
close all


t = 0 : 0.1 : 50;         % 1 kHz sample freq for 1 s
d = 0 : 5 : 50;           % 3 Hz repetition frequency
y = pulstran(t,d,'rectpuls',1);

plot(t,y)
xlabel 'Time (s)', ylabel Waveform