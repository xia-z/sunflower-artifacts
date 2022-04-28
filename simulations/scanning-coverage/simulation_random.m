% Created: 10/05/2021 by Qijia Shao
% Last modified: 10/05/2021 by Qijia Shao
% changes: random scan
    % 

close all; clc; clear;

% simulate the water surface using the model of sum of sines
num_sines = 3;
k = 2;
directions = [0; 3.8732; 5.9041];
amplitudes = [0.05; 0.000; 0.00];
wave_number = 2*3.14*[1/0.6; 1/0.1; 1/0.11]; %2pi/lambda
% % speeds V = lambda*f
incident_points = [0.1,0.1,0.1];
speeds = 2*3.14*[0.5;1.2;1.4];  %2pi*f
N = 100; % number of point for the wave
x = linspace(-0.2, 0.2, N); % range of the wave in the visualization
y = linspace(-0.1, 0.1, N);
[xx, yy] = meshgrid(x, y);

H = H_sum_of_sines(xx(:), yy(:), k, directions, amplitudes, wave_number, ...
    speeds, 0); % wave height

% global setting of the experiment
height = 2; % height of the laser
depth =  2; % depth of the underwater robot
refractive_index_air_water = 1.333;
laser_location = [0, 0, height]; 
beam_angle = 1.2; % out going beam angle of the laser
target_area_size = 2; % in meter 
robot_size = 0.0254; % one inch
searh_vector = -target_area_size/2: robot_size: target_area_size/2;
coverage_matrix = zeros(length(searh_vector), length(searh_vector));% representing if the cell is covered


spots = [-target_area_size/2,target_area_size/2];% underwater location (x,y)
pattern_points = []; % pattern locations at depth d


f1 = figure('name','Water Waves','NumberTitle','off');
set(0, 'DefaultLineLineWidth', 2);
s = surf(xx, yy, reshape(H, N, N)); % water surface
zlim([0, 0.05]);
axis equal

f2 = figure('name','Scanning pattern on the Surface','NumberTitle','off');
l2 = scatter(incident_points(:,1), incident_points(:,2)); % incident point location
hold on;
current2 = plot(0, 0, 'bo');
hold off;
xlim([-2, 2]);
ylim([-2, 2]);

f4 = figure('name','Scanning pattern under water','NumberTitle','off'); % show the trajectory
l4 = plot(spots(:, 1), spots(:, 2));  
hold on;
current4 = plot(0, 0, 'ro');
hold off;
xlim([-2, 2]);
ylim([-2, 2]);

f3 = figure('name','Scanning pattern under water','NumberTitle','off');
% l3 = plot(spots(:, 1), spots(:, 2), 'LineWidth', 5);  
l3 = scatter(spots(:, 1), spots(:, 2), 'filled', 'black'  );  % underwater spot location

size = (height+depth)* tan(beam_angle/360*pi); %Marker width in units of X underwater circle size

%Obtain the axes size (in axpos) in Points
currentunits = get(gca,'Units');
set(gca, 'Units', 'Points');
axpos = get(gca,'Position');
set(gca, 'Units', currentunits);
markerWidth = size/diff(xlim)*axpos(3); % Calculate Marker width in points
set(l3, 'SizeData', markerWidth^2); % set the marker size propotional to the axes
set(l2, 'SizeData', markerWidth^2);
hold on;
current3 = plot(0, 0, 'ro');
hold off;
xlim([-2, 2]);
ylim([-2, 2]);


% initialization for scanning
theta = 0;
t = 0;
r = 0;
n = 0; % number of points
m = 0; % number of points in the target area


grid_xy = -target_area_size/2:2*size:target_area_size/2;
[grid_x, grid_y] = meshgrid(grid_xy, grid_xy);
grid = [grid_x(:) grid_y(:)];
random_grid = randsample(length(grid),length(grid));
i = 0; % index indicate the number of samples

while i < length(grid)
    t = t + 0.05;
    i = i +1;
    pattern_point_x = grid(random_grid(i), 1);
    pattern_point_y = grid(random_grid(i), 2);
    pattern_point_z = -depth;
    pattern_points = [pattern_points; pattern_point_x, pattern_point_y];
    
   
    
    syms x y z;
    h_s = 2 * amplitudes(1) * ((((sin(wave_number(1)* x + speeds(1)*t) + 1) / 2) .^ k) - 1/2);
    [x0,y0,z0] = solve([x-laser_location(1) == (pattern_point_x - laser_location(1)) * (y - laser_location(2))/(pattern_point_y -laser_location(2)), ...
        x-laser_location(1) == (pattern_point_x - laser_location(1)) * (z - laser_location(3)) /(pattern_point_z -laser_location(3)),...
        z == h_s], [x,y,z], 'Real', true);

     incident_point_x = double([x0]);
     incident_point_y = double([y0]);
     incident_point_z = double([z0]);



    h = H_sum_of_sines(incident_point_x, incident_point_y, k, directions, amplitudes, wave_number, ...
    speeds, t); % wave height
    incident_points = [incident_points; incident_point_x, incident_point_y, incident_point_z];
  
    incident_direction = [incident_point_x, incident_point_y, incident_point_z] - laser_location;
    incident_direction = incident_direction/sqrt(sum(incident_direction.^2));
    
    H = H_sum_of_sines(xx(:), yy(:), k, directions, amplitudes, wave_number, speeds, t);  % whole water surface
    
    [Hx, Hy] = H_partial_sum_of_sines(incident_point_x, incident_point_y, k, directions, amplitudes, ...
        wave_number, speeds, t); % derivetives
    normal = [-Hx; -Hy; 1]; % normal at this point
    normal = normal / sqrt(sum(normal .^ 2));  % normalize the normal

 
%     calculate transmitted path
    incident_angle = acos(sum(-incident_direction' .* normal));
    sin_refracted_angle = sin(incident_angle) / refractive_index_air_water;
    refracted_angle = asin(sin_refracted_angle);  
    len = sin_refracted_angle / sin(incident_angle - refracted_angle);
    refracted = -normal + len * incident_direction';
    refracted = refracted/sqrt(sum(refracted.^2));
% %     refracted = incident_direction;
%     

      if ~isempty(refracted)
          a = incident_point_x + refracted(1) * (depth + h) / abs(refracted(3));
          b = incident_point_y + refracted(2) * (depth + h )/ abs(refracted(3));
        if refracted(3) >= 0 || a>1 || b >1
            spots = [spots; NaN, NaN];
        else
            spots = [spots;a,b];
            x_left = find(searh_vector <= a - size, 1, 'last');
            x_right = find(searh_vector <= a + size, 1, 'last');
            y_left = find(searh_vector <= b - size, 1, 'last');
            y_right = find(searh_vector <= b + size, 1, 'last');
            coverage_matrix(x_left:x_right, y_left: y_right) = 1;
        end
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
    
    l4.XData = spots(:, 1);
    l4.YData = spots(:, 2);
    current4.XData = spots(end, 1);
    current4.YData = spots(end, 2);
    
    
    if spots(end, 1)^2 + spots(end, 2)^2 <1
        m = m+1;
    end
%     M1(n) = getframe(f1);
%     M2(n) = getframe(f2);
%     M3(n) = getframe(f3);
    n = n + 1;
    
    
%     pause(0.0001);
end

coverage_count = nnz(coverage_matrix)
ratio = coverage_count/5929 % 5929 is when n2 = 1.0001
ratio_2 = coverage_count/4900 % 5929 when no wave but n2 = 1.33



