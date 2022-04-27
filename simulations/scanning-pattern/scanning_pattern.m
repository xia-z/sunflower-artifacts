% Simulation: modified archimedian spiral scan 

clear;
close all;

start = [20 20];
step = 0.5;
range = 90;
bucket_count = 5;
return_threshold = 1;
max_r = 10;

fig = figure;
hold on;
xlabel('x');
ylabel('y');
axis([-range range -range range]);
view(2);

theta = 0;
buckets = zeros(1, bucket_count);
while ishandle(fig)
	r = step * theta;
	x = r * cos(theta) + start(1);
	y = r * sin(theta) + start(2);
    delta = atan(step * (2 * pi) / r);
    color = 'red';
	if (sqrt((r * cos(theta))^2 + (r * sin(theta))^2)) >= max_r 
		break
	end
    if sqrt(x^2 + y^2) >= range
        theta_c = 2 * atan(y/(x + sqrt(x^2 + y^2)));
        x = range * cos(theta_c);
		y = range * sin(theta_c);
        delta = atan(step * (2 * pi) / sqrt((x - start(1))^2 + (y - start(2))^2));
        for i = 1:bucket_count
            minRange = (2*pi)/bucket_count * (i - 1);
            maxRange = (2*pi)/bucket_count * (i);
            angle = theta_c + pi;
            if angle > minRange && angle <= maxRange
                buckets(i) = 1;
            end
        end
        color = 'green';
    end
    theta = theta + delta;
    b = viscircles([x y], step * pi, 'Color', color);
    drawnow;
    if length(buckets(buckets == 0)) <= 1
		break
    end
end
