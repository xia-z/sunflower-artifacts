close all; clc; clear;
% simulate the water surface using the model of sum of sines
num_sines = 3;
k = 2;

directions = [0; 3.8732; 5.9041];
amplitudes = [0.05; 0.00; 0.00];
wave_number = 2*3.14159*[1/0.6; 1/0.1; 1/0.11]; %2pi/lambda
% % speeds V = lambda*f
incident_points = [0,0,0];
speeds = 2*3.14*[1;1.2;1.4];  %2pi*f



N = 100;

height = 2;
depth =  0.1;
refractive_index_air_water = 1.33;
laser_location = [0, 0, height]; % spiral scan 

x = linspace(-0.2, 0.2, N);
y = linspace(-0.1, 0.1, N);
[xx, yy] = meshgrid(x, y);

H = H_sum_of_sines(xx(:), yy(:), k, directions, amplitudes, wave_number, ...
    speeds, 0);
spots = [0, 0];

f1 = figure('name','Water Waves','NumberTitle','off');
set(0, 'DefaultLineLineWidth', 2);
s = surf(xx, yy, reshape(H, N, N)); % water surface
zlim([0, 0.4]);
axis equal

f2 = figure('name','Scanning pattern on the Surface','NumberTitle','off');
l2 = plot(incident_points(:,1), incident_points(:,2)); % incident point location
hold on;
current2 = plot(0, 0, 'bo');
hold off;
xlim([-2, 2]);
ylim([-2, 2]);

f3 = figure('name','Scanning pattern under water','NumberTitle','off');
l3 = plot(spots(:, 1), spots(:, 2));
hold on;
current3 = plot(0, 0, 'ro');
hold off;
xlim([-2, 2]);
ylim([-2, 2]);

n = 1;


% incident_angles = [0];
% normals = [];

for t = 0:0.05:10
    
    % spiral scan
    r = 1/10*t; % radius
    incident_point_x = r* cos(2*3.14*t);
    incident_point_y = r* sin(2*3.14*t);
    
    h = H_sum_of_sines(incident_point_x, incident_point_y, k, directions, amplitudes, wave_number, ...
    speeds, t); % wave height
    incident_point_z = h;
%     incident_point_x = 0;
%     incident_point_y = 0;
    incident_points = [incident_points; incident_point_x, incident_point_y, 0]; % just for compute the incident direction
    
    incident_direction = [incident_point_x, incident_point_y, incident_point_z] - laser_location;
    incident_direction = incident_direction/sqrt(sum(incident_direction.^2));
    
    H = H_sum_of_sines(xx(:), yy(:), k, directions, amplitudes, wave_number, speeds, t);  % whole water surface
    
    [Hx, Hy] = H_partial_sum_of_sines(incident_point_x, incident_point_y, k, directions, amplitudes, ...
        wave_number, speeds, t); % derivetives
    normal = [-Hx; -Hy; 1]; % normal at this point
    normal = normal / sqrt(sum(normal .^ 2));  % normalize the normal

 
    % calculate transmitted path
    incident_angle = acos(sum(-incident_direction' .* normal));
    sin_refracted_angle = sin(incident_angle) / refractive_index_air_water;
    refracted_angle = asin(sin_refracted_angle);  
    len = sin_refracted_angle / sin(incident_angle - refracted_angle);
    refracted = -normal + len * incident_direction';
    refracted = refracted/sqrt(sum(refracted.^2));
%     refracted = incident_direction;
    
    if refracted(3) >= 0
        spots = [spots; NaN, NaN];
    else
        spots = [spots; ...
            incident_point_x + refracted(1) * (depth ) / abs(refracted(3)), ...
            incident_point_y + refracted(2) * (depth ) / abs(refracted(3))];
    end
    

    s.ZData = reshape(H, N, N);
    
    l2.XData =  incident_points(:,1);
    l2.YData =  incident_points(:,2);
    current2.XData = incident_points(end, 1);
    current2.YData = incident_points(end, 2);
    
    l3.XData = spots(:, 1);
    l3.YData = spots(:, 2);
    current3.XData = spots(end, 1);
    current3.YData = spots(end, 2);
    
    M1(n) = getframe(f1);
    M2(n) = getframe(f2);
    M3(n) = getframe(f3);
    n = n + 1;
    
    
    pause(0.05)
end

% outputVideo1 = VideoWriter('wave');
% outputVideo1.FrameRate = 20;
% open(outputVideo1);
% for i = 1:size(M1, 2)
%     writeVideo(outputVideo1, M1(i));
% end
% close(outputVideo1);

