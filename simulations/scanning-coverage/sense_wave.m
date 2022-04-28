function h = sense_wave(x, y, t, wave_parameters, sensing_parameters)
num_sensors = size(x, 1);
h = zeros(num_sensors, 1);
k = wave_parameters.k;
D = wave_parameters.D;
A = wave_parameters.A;
omega = wave_parameters.omega;
phi = wave_parameters.phi;
distance_tx_water = wave_parameters.distance_tx_water;

sensing_delay = sensing_parameters.sensing_delay;
sensing_noise_sigma = sensing_parameters.sensing_noise_sigma;

for i = 1:num_sensors
    h(i) = H_sum_of_sines(x(i), y(i), k, D, A, omega, phi, t + sensing_delay * (i-1));
end
h = h - distance_tx_water;
h = h + normrnd(0, sensing_noise_sigma, num_sensors, 1);
end