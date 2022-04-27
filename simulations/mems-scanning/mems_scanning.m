% Simulation: MEMS mirror scanning trajectories

% Reset
close all;
clearvars;

% Scanning vars
mechanical_max = 3; % max mechanical angle, full mechanical FOV is ±mechanical_max, deg
lens_diameter = 10; % scanning diameter, mm
beam_diameter = 1; % mm
start_xy = [2 0]; % starting points [x y]
bucket_threshold = 5; % Number of missed buckets allowed
%trajectory_resolution = 0.001; % interpolation resolution between points
scale = 500; % just for display purposes
show_scan = true; % should show display
return_threshold = 0.1; % how close to [0 0] does the return trajectory need to get
scan_frequency = 100; % hz

% Scanning calculations
theta = 0; 
target_found = false;
scanning_finished = false;
spacing = beam_diameter/(2 * pi); % Spiral spacing
bucket_delta = atan(beam_diameter/(lens_diameter/2));
buckets = zeros(1, int8(2*pi/bucket_delta)); % Edge buckets
distance = (lens_diameter/2)/tan(mechanical_max * 2 * pi/180); % distance from mirror to scanning area
y_angles = atan(start_xy(2)/distance) * 180/pi * 0.5;
x_angles = atan(start_xy(1)/distance) * 180/pi * 0.5;
    
% Setup figures
fig = figure('Position', [0 0 3 * scale 1 * scale], 'visible', 'off');
movegui(fig, 'center');
if show_scan
	
    % Setup first plot
    subplot(1, 3, 1);
    hold on;
    grid on;
    title(sprintf('Lens Hitting Position (%.2f mm from MEMS)', distance));
    xlabel('X (mm)');
    ylabel('Y (mm)');
    pbaspect([1 1 1]);
    plot(lens_diameter/2 * cos(0:pi/50:2*pi), lens_diameter/2 * sin(0:pi/50:2*pi), '-', 'Color', 'blue'); % Steering boundary
    axis([-lens_diameter/2 lens_diameter/2 -lens_diameter/2 lens_diameter/2]); % Set axis
    c_x = 0;
    c_y = 0;
    c = plot(c_x, c_y, '-', 'Color', 'black');
    hold off;

    % Setup second plot
    subplot(1, 3, 2);
    hold on;
    grid on;
    title('MEMS X-Axis Mechanical Angle');
    xlabel('Step');
    ylabel('Angle (°)');
    pbaspect([1 1 1]);
    ylim([-mechanical_max mechanical_max]); % Set axis
    m_px = plot(0, x_angles, '-');
    hold off;

    % Setup third plot
    subplot(1, 3, 3);
    hold on;
    grid on;
    title('MEMS Y-Axis Mechanical Angle');
    xlabel('Step');
    ylabel('Angle (°)');
    pbaspect([1 1 1]);
    ylim([-mechanical_max mechanical_max]); % Set axis
    m_py = plot(0, y_angles, '-');
    hold off;

    % Draw
    drawnow;
    set(fig, 'visible', 'on');
end

% Start the scan
i = 0;
while ~target_found && ~scanning_finished && ishandle(fig)
	
	% Calculate radius for spiral
	r = spacing * theta;
	
	% Compute cartesian coordiantes of beam's center
	x = r * cos(theta) + start_xy(1);
	y = r * sin(theta) + start_xy(2);
	
	% Beam line color
	line_color = 'red';
    
	% Determine if the coordaintes are outside the boundary
	if sqrt(x^2 + y^2) >= lens_diameter/2
		
		% Calculate theta relative to the center
		center_theta = 2 * atan(y/(x + sqrt(x^2 + y^2)));
		
		% Assign each boundary angle to a bucket
		bucket = floor((center_theta + pi)/bucket_delta) + 1;
		buckets(bucket) = 1;
		
		% Calculate actual x and y coords
		x = lens_diameter/2 * cos(center_theta);
		y = lens_diameter/2 * sin(center_theta);

		% Set different line color
		line_color = 'black';
		
		% Update theta on the boundaries
		%r_center = lens_diameter/2;
		%r_relative = r;
		%x_offset = start(1);
		%y_offset = start(2);
		%delta_theta_center = 2 * acos((r_center^2 + x_offset^2 + y_offset^2 - r_relative^2)/(2 * r_center * sqrt(x_offset^2 + y_offset^2)));
		%delta_theta_relative = 2 * asin(r_center * sin(delta_theta_center / 2)));
		theta = theta + atan(beam_diameter / sqrt((x - start_xy(1))^2 + (y - start_xy(2))^2));
		%theta = theta + atan(beam_diameter / r);
	else
		
		% Update theta
		theta = theta + atan(beam_diameter / r);
    end
    
	% Save mirror scanning angles
	x_angles = [x_angles atan(x/distance) * 180/pi * 0.5];
	y_angles = [y_angles atan(y/distance) * 180/pi * 0.5];
	
	
    % Draw the beam
    if show_scan
        subplot(1,3,1);
        b = viscircles([x y], beam_diameter / 2, 'LineWidth', 1, 'LineStyle', '-', 'Color', line_color);
  
        % Set center data
        c.XData = [0 x];
        c.YData = [0 y];
	
		% Set x,y data
        m_px.XData = 1:length(x_angles);
        m_px.YData = x_angles;
        m_py.XData = 1:length(y_angles);
        m_py.YData = y_angles;
	
        % Draw
        drawnow;
    end
	i = i + 1;
	% End if the number of missed buckets is below the threshold
	if length(buckets(buckets == 0)) <= bucket_threshold
		scanning_finished = true;
    end
end


% Return to zero
last_r = r;
returned = false;
while ~returned
	r = 2 * last_r - spacing * theta;
	x = r * cos(theta);
	y = r * sin(theta);
	if sqrt(x^2 + y^2) >= lens_diameter/2
		center_theta = 2 * atan(y/(x + sqrt(x^2 + y^2)));
		x = lens_diameter/2 * cos(center_theta);
		y = lens_diameter/2 * sin(center_theta);
	end
	if abs(x) <= return_threshold && abs(y) <= return_threshold
		returned = true;
	end
	x_angles = [x_angles atan(x/distance) * 180/pi * 0.5];
	y_angles = [y_angles atan(y/distance) * 180/pi * 0.5];
	theta = theta + atan(beam_diameter / r);
	
	% Draw the beam
	if show_scan
        subplot(1,3,1);
        b = viscircles([x y], beam_diameter / 2, 'LineWidth', 1, 'LineStyle', '-', 'Color', 'yellow');
  
        % Set center data
        c.XData = [0 x];
        c.YData = [0 y];
	
		% Set x,y data
        m_px.XData = 1:length(x_angles);
        m_px.YData = x_angles;
        m_py.XData = 1:length(y_angles);
        m_py.YData = y_angles;
	
        % Draw
        drawnow;
	end
	i = i + 1;
end

i
fig3 = figure;
plot(x_angles, y_angles);

% Convert to time
x_seconds = 0:1/scan_frequency:length(x_angles)/scan_frequency - 1/scan_frequency;
y_seconds = 0:1/scan_frequency:length(y_angles)/scan_frequency - 1/scan_frequency;


% New figure for interpolated trajectories
fig2 = figure('Position', [0 0 2 * scale 1 * scale]);
subplot(1, 2, 1);
hold on;
grid on;
title('MEMS X-Axis Trajectory');
xlabel('Step');
ylabel('Angle (°)');
pbaspect([1 1 1]);
ylim([-mechanical_max mechanical_max]); % Set axis
x_interp = x_angles;%interp1(1:length(x_angles), x_angles, 1:trajectory_resolution:length(x_angles), 'spline');
plot(x_interp);
hold off;
subplot(1, 2, 2);
hold on;
grid on;
title('MEMS Y-Axis Trajectory');
xlabel('Step');
ylabel('Angle (°)');
pbaspect([1 1 1]);
ylim([-mechanical_max mechanical_max]); % Set axis
y_interp = y_angles;%interp1(1:length(y_angles), y_angles, 1:trajectory_resolution:length(y_angles), 'spline');
plot(y_interp);
hold off;


%{
% Final trajectories
fig2 = figure('Position', [0 0 2 * scale 1 * scale]);
subplot(1, 2, 1);
hold on;
grid on;
title('MEMS X-Axis Trajectory');
xlabel('Time (s)');
ylabel('Angle (°)');
pbaspect([1 1 1]);
ylim([-mechanical_max mechanical_max]); % Set axis
plot(x_seconds, x_angles);
hold off;
subplot(1, 2, 2);
hold on;
grid on;
title('MEMS Y-Axis Trajectory');
xlabel('Time (s)');
ylabel('Angle (°)');
pbaspect([1 1 1]);
ylim([-mechanical_max mechanical_max]); % Set axis
plot(y_seconds, y_angles);
hold off;
%}