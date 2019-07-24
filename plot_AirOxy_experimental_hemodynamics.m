% Plot the experimental data for either air or oxy cases with data averaged over region, animal and trial number

figure;
hold all
plot(tim, HbO_oxy);
plot(tim, HbT_oxy);
plot(tim, HbR_oxy);
legend('HbO','HbT','HbR');
title('Oxy')

figure;
hold all
plot(tim, HbO_air);
plot(tim, HbT_air);
plot(tim, HbR_air);
legend('HbO','HbT','HbR');
title('Air')