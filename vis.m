plot(0 : length(energy) - 1, energy, 'b')
hold on;
plot(0 : length(Elja_energy) - 1, Elja_energy, 'g')
xlabel('iteration')
ylabel('energy')
title('Energy')
legend('\alpha\beta-swap', '\alpha-expansion');%, 'Location', 'NorthWest')
print('energy_2Elja', '-depsc2', '-r300');
    
plot(0 : 2, [0, time(1:2)], 'r')
hold on;
plot(0 : length(Elja_time), [0, Elja_time], 'g')
xlabel('iteration')
ylabel('time, sec')
title('Time')
legend('\alpha\beta-swap', '\alpha-expansion', 'Location', 'NorthWest')
print('time_2Elja', '-depsc2', '-r300');

energy = energy(energy > 0);
plot(0 : length(energy) - 1, energy, 'b')
xlabel('iteration')
ylabel('energy')
legend('energy')
print('energy_mC', '-depsc2', '-r300');

plot(0 : 5, [0, time(1 : 5)], 'r')
xlabel('iteration')
ylabel('time')
legend('time', 'Location', 'NorthWest')
print('time_mC', '-depsc2', '-r300');