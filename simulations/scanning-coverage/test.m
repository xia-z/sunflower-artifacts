% simulate the water surface using the model of sum of sines
num_sines = 3;
k = 2;

% directions = [0; 1/30 * pi];
% amplitudes = [0.03; 0.02];
% frequencies = [2 * pi / 0.1; 2 * pi / 0.1];
% speeds = [10 / 7.5 * 2 * pi; 10 / 7.5 * 2 * pi];
% directions = 2 * pi * rand(num_sines, 1)
% directions = zeros(num_sines, 1) + [0; 2 * pi * rand(num_sines - 1, 1)]

% amplitudes = [0.01; abs(randn(num_sines - 1, 1) * 0.002)];
% amplitudes = rand(num_sines, 1) * 0.05
% frequencies = rand(num_sines, 1) * 10
% frequencies = [2 * pi / 0.2; rand(num_sines - 1, 1) * 100]
% speeds = rand(num_sines, 1)
% speeds = 10 / 7.5 * 2 * pi * (ones(num_sines, 1) + [0; randn(num_sines - 1, 1)])

directions = [0; 3.8732; 5.9041];
amplitudes = [0.0100; 0.000; 0.000];
frequencies = [31.4159; 98.4349; 94.5579];
% speeds = [8.3776; 12.4452; 16.9731];
speeds = [8.3776; 10.4452; 12.9731];
N = 100;

H_Rx = 0.5;
depth = 0.2;
refractive_index_air_water = 1.33;

x = linspace(-0.2, 0.2, N);
y = linspace(-0.1, 0.1, N);
[xx, yy] = meshgrid(x, y);

H = H_sum_of_sines(xx(:), yy(:), k, directions, amplitudes, frequencies, ...
    speeds, 0);
trajectory = [0, 0];
spots = [0, 0];

f1 = figure(1);
set(0, 'DefaultLineLineWidth', 2);
% subplot(1, 2, 1);
s = surf(xx, yy, reshape(H, N, N));
zlim([0, 0.05]);
% sfh1.Position = [0, 0, 100, 100];
axis equal
f2 = figure(2);
% subplot(1, 2, 2);
% set(sfh2, 'Position', [100, 0, 100, 100]);
l = plot(trajectory(:, 1), trajectory(:, 2));
hold on;
current = plot(0, 0, 'ro');
hold off;
axis equal
xlim([-0.5, 0.5]);
ylim([-0.5, 0.5]);

f3 = figure(3);
l3 = plot(spots(:, 1), spots(:, 2));
hold on;
current3 = plot(0, 0, 'ro');
hold off;
xlim([-0.1, 0.1]);
ylim([-0.1, 0.1]);

n = 1;

% incident_angles = [0];
% normals = [];

for t = 0:0.05:10
    H = H_sum_of_sines(xx(:), yy(:), k, directions, amplitudes, frequencies, ...
    speeds, t);
    
    [Hx, Hy] = H_partial_sum_of_sines(0, 0, k, directions, amplitudes, ...
        frequencies, speeds, t);
    normal = [-Hx; -Hy; 1];
    normal = normal / sqrt(sum(normal .^ 2));  % normalize the normal
%     normals = [normals; normal'];
    len = sqrt(2 + 2 * (2 * normal(3)^2 - 1));
    reflected = len * normal - [0; 0; 1];
    h = H_sum_of_sines(0, 0, k, directions, amplitudes, frequencies, ...
        speeds, t);
    
    if reflected(3) <= 0
        trajectory = [trajectory; NaN, NaN];
    else
        trajectory = [trajectory; ...
            reflected(1) * (H_Rx - h) / reflected(3), ...
            reflected(2) * (H_Rx - h) / reflected(3)];
    end
    
    % calculate transmitted path
    incident_angle = acos(sum([0; 0; 1] .* normal));
%     incident_angles = [incident_angles; incident_angle];
    sin_refracted_angle = sin(incident_angle) / refractive_index_air_water;
    refracted_angle = asin(sin_refracted_angle);
    len = sin_refracted_angle / sin(incident_angle - refracted_angle);
    refracted = -normal + len * [0; 0; -1];
    if refracted(3) >= 0
        spots = [spots; NaN, NaN];
    else
        spots = [spots; ...
            refracted(1) * (depth + h) / abs(refracted(3)), ...
            refracted(2) * (depth + h) / abs(refracted(3))];
    end
    
    % calculate transmitted path from water to air
%     incident_angle = acos(sum([0; 0; 1] .* normal));
%     sin_refracted_angle = sin(incident_angle) * refractive_index_air_water;
%     refracted_angle = asin(sin_refracted_angle);
%     len = sin(pi - refracted_angle) / sin(refracted_angle - incident_angle);
%     refracted = len * [0; 0; 1] - normal;
%     if refracted(3) <= 0
%         spots = [spots; NaN, NaN];
%     else
%         spots = [spots; ...
%             refracted(1) * (depth - h) / refracted(3), ...
%             refracted(2) * (depth - h) / refracted(3)];
%     end
    
    s.ZData = reshape(H, N, N);
    l.XData = trajectory(max(1, end - 100):end, 1);
    l.YData = trajectory(max(1, end - 100):end, 2);
    current.XData = trajectory(end, 1);
    current.YData = trajectory(end, 2);
    l3.XData = spots(:, 1);
    l3.YData = spots(:, 2);
    current3.XData = spots(end, 1);
    current3.YData = spots(end, 2);
    
    M1(n) = getframe(f1);
    M2(n) = getframe(f2);
    M3(n) = getframe(f3);
    n = n + 1;
    
%     subplot(1, 2, 2);
%     l = plot(trajectory(:, 1), trajectory(:, 2));
%     axis equal
%     xlim([-10, 10]);
%     ylim([-10, 10]);
    
    pause(0.05)
end

% outputVideo1 = VideoWriter('wave');
% outputVideo1.FrameRate = 20;
% open(outputVideo1);
% for i = 1:size(M1, 2)
%     writeVideo(outputVideo1, M1(i));
% end
% close(outputVideo1);
% 
% outputVideo2 = VideoWriter('trace');
% outputVideo2.FrameRate = 20;
% open(outputVideo2);
% for i = 1:size(M2, 2)
%     writeVideo(outputVideo2, M2(i));
% end
% close(outputVideo2);
% 
% outputVideo3 = VideoWriter('refracted');
% outputVideo3.FrameRate = 20;
% open(outputVideo3);
% for i = 1:size(M3, 2)
%     writeVideo(outputVideo3, M3(i));
% end
% close(outputVideo3);
