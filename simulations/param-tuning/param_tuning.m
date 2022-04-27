% Simulation: Parameter tuning for optimizing sensing distance

close all;
clc;

d_tx = 2/1000; % diameter of TX aperture [m]
d_rr = 25.4/1000; % diameter of retoreflector [m]
d_rx = 40/1000; % diameter of fiber ring [m]
drange = 0:1:200; % full distances to consider [m]
trange = [1 2 4 6 8 10 20 30]; % half-angle, [mrads]
P_tx = 80; % tx power [mW]
P_noise = 1; % noise level, experimentally measured [mW]

y_max = 0;
for t = trange/1000
    x = [];
    y = [];
    for d = drange
        P_rr = P_tx * min(d_rr^2/(d_tx+tan(t)*d)^2, 1);
        P_rx = P_rr * min(d_rx^2/(d_rr+tan(t)*d)^2, 1);
        SNR = 10*log10(P_rx/P_noise);
        x = [x d];
        y = [y SNR];
        if y(end) > y_max
            y_max = y(end);
        end
    end
    hold on;
    plot(x, y, 'DisplayName', sprintf('%.1f°', 2*t*180/pi));
    hold off;
end

title(sprintf('P_0 = %i mW', P_tx))
lgd = legend;
lgd.Title.String = 'Div (Full-Angle)';
ylabel('Sensing SNR (dB)');
xlabel('Drone-to-AUV distance (m)');
% set(gca, 'YScale', 'log');
% set(gca, 'XScale', 'log');
line = yline(3, 'Label', '3dB');
line.Annotation.LegendInformation.IconDisplayStyle = 'off';
ylim([-3 y_max]);
