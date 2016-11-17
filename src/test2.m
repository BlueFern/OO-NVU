
factor=[2 1 0.5 0.1 0.01 ];
topamplitude=[14.25 15.54 16.06 16.36 16.42 ];
bottomamplitude=[13.89 14.72 14.77 14.82 14.83 ];
frequency=[12/90 6/65 7/90 7/100 6/87 ];

figure(1)
subplot(1,2,1)
scatter(factor, topamplitude-bottomamplitude, 'fill')
h = lsline;
lsline
set(h,'linewidth',2, 'color', 'k')
title('amplitude [\mum]')
xlabel('G_C_a * factor [-]')

subplot(1,2,2)
scatter(factor, frequency, 'fill')
h = lsline;
lsline
set(h,'linewidth',2, 'color', 'k')
title('frequency [Hz]')
xlabel('G_C_a * factor [-]')
