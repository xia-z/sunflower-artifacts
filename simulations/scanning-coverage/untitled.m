% x = -0.1:0.00254:0.1;
% y = -0.1:0.00254:0.1;
% 
% point = [0.4, 0.3];
% x_location_left = find(x <= point(1)-0.04, 1, 'last');
% x_location_right = find(x <= point(1)+0.04, 1, 'last');
% 
% y_location_left = find(y <= point(1)-0.04, 1, 'last');
% y_location_right = find(y <= point(1)+0.04, 1, 'last');
% 
% coverage_matrix = zeros(length(x), length(y));
% coverage_matrix(x_location_left:x_location_right, y_location_left:y_location_right) = 1;
% 
% x(x_location)
% y(y_location)
% N = nnz(coverage_matrix) 

% define some grids
step_size = 0.05;
target_area_size = 2;
points = [-target_area_size/2,target_area_size/2];
x = -target_area_size/2:step_size:target_area_size/2;
y = -target_area_size/2:step_size:target_area_size/2;
[xx, yy] = meshgrid(x, y);
point1 = [xx(:) yy(:)];
i = 0;

f4 = figure('name','Scanning pattern under water','NumberTitle','off'); % show the trajectory
l4 = scatter(points(:,1),points(:,2));  
xlim([-2, 2]);
ylim([-2, 2]);


mm = randsample(length(point1),length(point1));
