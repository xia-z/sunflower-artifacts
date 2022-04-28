close all; clc; clear;
% simulate the water surface using the model of sum of sines
num_sines = 3;
k = 2;

directions = [0; 3.8732; 5.9041];
amplitudes = [0.05; 0.03; 0.02];
wave_number = 2*3.14*[1/0.6; 1/0.1; 1/0.11]; %2pi/lambda
% % speeds V = lambda*f
incident_points = [0,0,0];
speeds = 2*3.14*[1;1.2;1.4];  %2pi*f

N = 100;

height = 2;
depth =  0.1;
refractive_index_air_water = 1.33;
laser_location = [0, 0, height]; % spiral scan 
beam_angle = 3;  % 1.2

x = linspace(-0.2, 0.2, N);
y = linspace(-0.1, 0.1, N);
[xx, yy] = meshgrid(x, y);

H = H_sum_of_sines(xx(:), yy(:), k, directions, amplitudes, wave_number, ...
    speeds, 0);
spots = [0, 0];
pattern_points = [0,0];

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

f4 = figure('name','Scanning pattern under water','NumberTitle','off');
l4 = plot(spots(:, 1), spots(:, 2));  
hold on;
current4 = plot(0, 0, 'ro');
hold off;
xlim([-2, 2]);
ylim([-2, 2]);

f3 = figure('name','Scanning pattern under water','NumberTitle','off');
% l3 = plot(spots(:, 1), spots(:, 2), 'LineWidth', 5);  
% l3 = viscircles([spots(:,1), spots(:,2)], 1) ;
l3 = scatter(spots(:, 1), spots(:, 2));  % underwater spot location

size = (height+depth)* tan(beam_angle/360*pi); %Marker width in units of X underwater circle size
% size = 0.5;
%Obtain the axes size (in axpos) in Points
currentunits = get(gca,'Units');
set(gca, 'Units', 'Points');
axpos = get(gca,'Position');
set(gca, 'Units', currentunits);
markerWidth = size/diff(xlim)*axpos(3); % Calculate Marker width in points
set(l3, 'SizeData', markerWidth^2);
set(l2, 'SizeData', markerWidth^2);
hold on;
current3 = plot(0, 0, 'ro');
hold off;
xlim([-2, 2]);
ylim([-2, 2]);




% initialization for scanning
theta = 0;
step = size/pi;
t = 0;
r = 0;
n = 0; % number of points
m = 0; % number of points in the target area


while(r <= 1)
    t = t + 0.05;
    
    % spiral scan
%     r = 1/10*t; % radius
    r = step * theta;
	pattern_point_x = r * cos(theta) ;
	pattern_point_y = r * sin(theta) ;
    delta = atan(step * (2 * pi) / r);
    theta = theta + delta;
%     pattern_point_x = r* cos(2*3.14*t);
%     pattern_point_y = r* sin(2*3.14*t);
%     pattern_point_x = 0.5;
%     pattern_point_y = 0;
    pattern_point_z = -depth;
    pattern_points = [pattern_points; pattern_point_x, pattern_point_y];
    
    h_s = 2 * amplitudes(1) * ((((sin(wave_number(1)* x + speeds(1)*t) + 1) / 2) .^ k) - 1/2);
    
    syms x y z;
%     equals = [x-laser_location(1) == (pattern_point_x - laser_location(1)) * (y - laser_location(2))/(pattern_point_y -laser_location(2)),...
%         x-laser_location(1) == (pattern_point_x - laser_location(1)) * (z - laser_location(3)) /(pattern_point_z -laser_location(3)),...
%         z == h_s]; 
    
%     e1 = x-laser_location(1) -(patten_point_x - laser_location(1)) * (y - laser_location(2))/(patten_point_y -laser_location(2));
%     e2 = x-laser_location(1) - (patten_point_x - laser_location(1)) * (z - laser_location(3)) /(patten_point_z -laser_location(3));
%     e3 = z - H_sum_of_sines(x, y, k, directions, amplitudes, wave_number, speeds, t);

    [x0,y0,z0] = solve([x-laser_location(1) == (pattern_point_x - laser_location(1)) * (y - laser_location(2))/(pattern_point_y -laser_location(2)), ...
        x-laser_location(1) == (pattern_point_x - laser_location(1)) * (z - laser_location(3)) /(pattern_point_z -laser_location(3)),...
        z == h_s], [x,y,z], 'Real', true);

     incident_point_x = double([x0]);
     incident_point_y = double([y0]);
     incident_point_z = double([z0]);



    h = H_sum_of_sines(incident_point_x, incident_point_y, k, directions, amplitudes, wave_number, ...
    speeds, t); % wave height
%     incident_point_z = h;

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
        if refracted(3) >= 0 || (a^2+b^2)>1
            spots = [spots; NaN, NaN];
        else
            spots = [spots;a,b];
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
    
    
    pause(0.001)
end

set(gca, 'Visible', 'off')
saveas(f3,'test.png')
Image  = imread('test.png');
BW     = im2bw(Image);
nBlack = 1021125 - sum(BW(:));
nBlack/110611


% 
% outputVideo1 = VideoWriter('wave');
% outputVideo1.FrameRate = 20;
% open(outputVideo1);
% for i = 1:size(M1, 2)
%     writeVideo(outputVideo1, M1(i));
% end
% close(outputVideo1);

