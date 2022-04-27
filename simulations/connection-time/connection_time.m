% Simulation: Computes percentage of hit/miss given a diverging laser beam
% and water wave

clearvars;
clc;
close all;

% Params
syms a b c d;
p = [];

% Beam characteristics
p.div = 1.2; % degrees
p.offset = 0+1e-12;

% Air
p.h = 4;
p.n1 = 1;

% Water
p.n2 = 1.33;
p.d = 4;
p.rx = 0.0254;

% Wave
p.A = 1.5; % peak to zero
p.l = 33.8; % wavelength
p.T = 5.7; % 1/f

% Equations
p.Wave = d*sin(a/p.l - b/p.T);
p.TanSlope = diff(p.Wave, a);
p.Intercept = p.h * tan(c) * (p.h - p.Wave)/p.h == a;

% Plotting
p.ti = 3;
p.ni = 1;
p.dt = 0.05;
p.dx = 0.05;
p.x = -1*pi:p.dx:1*pi;
p.t = 0:p.dt:4*p.T*pi;
p.plot = true;
p.show_norm = false;
p.show_tan = false;
p.show_origin = false;
p.show_hypotenuse = false;

% Calculations
hits = [];

% Plot for one full loop
for t = p.t
	
	% Trace left and right rays
	lpts = traceray(p, t, -p.div/2+p.offset, p.A);
	rpts = traceray(p, t, p.div/2+p.offset, p.A);

	% Trace center point of beam
	cpts = traceray(p, t, p.offset, 0.0001);
	rx_c = cpts.final(1);
	rx_l = rx_c-p.rx/2;
	rx_r = rx_c+p.rx/2;

	% Determine hit/miss
	if rpts.final(1) >= rx_l && lpts.final(1) <= rx_r
		hits = [hits 1];
	else
		hits = [hits 0];
	end

	% Plot
	if p.plot

		% Set up figure
		clf(gcf);
		hold on;
		axis equal;
		grid on;
		xlim([p.x(1), p.x(end)]);
		ylim([-p.d, p.h]);

		% Plot wave
		y_wave = subs(p.Wave, {a, b, d}, {p.x, t, p.A});
		plot(p.x, y_wave, 'color', 'blue');

		% Plot beams
		plotray(p, t, lpts);
		plotray(p, t, rpts);

		% Plot RX
		scatter(rx_l, -p.d, 100, 'black', 'marker', '|');
		scatter(rx_r, -p.d, 100, 'black', 'marker', '|');

		% Plot hit
		if hits(end) == 1
			title('Hit');
		else
			title('Miss');
		end

		% Draw
		hold off;
		drawnow;
	end
end
hold off;

% Compute hit duration
ts = [0];
% figure;
for i = 2:length(p.t)
	if hits(i) == 1 && hits(i-1) == 0
% 		xline(p.t(i), label='Start', LabelHorizontalAlignment='left');
		ts = [ts p.t(i)];
	elseif hits(i) == 0 && hits(i-1) == 1
% 		xline(p.t(i-1), label='Stop');
		ts = [ts p.t(i-1)];
	end
end
ts = [ts p.t(end)];
dts = diff(ts);
% dts = dts(2:end-1);
states = abs(rem(1:length(dts), 2)-1*~hits(1));
% states = states(2:end-1);
on = dts(states == 1);
off = dts(states == 0);
fprintf('HIT min=%.2fs max=%.2fs mean=%.2fs', min(on), max(on), mean(on));
fprintf('MISS min=%.2fs max=%.2fs mean=%.2f ', min(off), max(off), mean(off));
fprintf('HIT=%.2f%%, MISS=%.2f%%, t=%.2fs', sum(on)/sum([on off])*100, sum(off)/sum([on off])*100, p.t(end))


function points = traceray(p, t, alpha, A)
	syms a b c d;
	
	% Beam starting point
	x_0 = 0;
	y_0 = p.h;
	
	% Steering angle in rads
	th_steering = alpha*pi/180; 

	% x/y intercept of beam with wave
	x_w = vpasolve(subs(p.Intercept, {b, c, d}, {t, th_steering, A}), a);
	y_w = subs(p.Wave, {a, b, d}, {x_w, t, A});

	% Slopes at intercept relative to origin's reference frame
	m_tangent = subs(p.TanSlope, {a, b, d}, {x_w, t, A}); 
	m_normal = -1/m_tangent; 
	m_incident = -1/tan(th_steering); 
	
	% Incident/refracted angles of beam at intercept relative to normal
	th_incident = (atan((m_normal-m_incident)/(1+m_normal*m_incident))) * -m_normal/abs(m_normal);
	th_refracted = asin(sin(th_incident)*p.n1/p.n2);
	
	% Normal's rotation relative to origin's reference frame
	th_normal_rotation = abs(atan(1/m_normal));
	
	% Refracted beam's angle at intercept relative to origin's frame
	th_refracted_origin = (th_normal_rotation - th_refracted);
	if m_normal > 0 % I don't know why this works, but it works...
		th_refracted_origin = (th_refracted - th_normal_rotation);
	end

	% Beam end point
	x_1 = x_w + (p.d+y_w)*tan(th_refracted_origin);
	y_1 = -p.d;

	% Return critical points
	points = [];
	points.initial = [x_0, y_0];
	points.intercept = [x_w, y_w];
	points.final = [x_1, y_1];
end

function plotray(p, t, pts)

	% Extract points
	initial = pts.initial;
	intercept = pts.intercept;
	final = pts.final;
	
	% Background lines
	if p.show_origin
		xline(double(intercept(1)));
		yline(double(intercept(2)));
	end
	if p.show_hypotenuse
		plot([initial(1), final(1)], [initial(2), final(2)], 'color', 'green')
	end

	% Draw origin, intercept, and final points
	scatter(initial(1), initial(2), 100, 'black', 'marker', '*')
	scatter(intercept(1), intercept(2), 500, 'black', 'marker', '.');
	scatter(final(1), final(2), 100, 'black', 'marker', 'x');

	% Draw rays
	plot([initial(1), intercept(1)], [initial(2), intercept(2)], 'color', 'red'); 
	plot([intercept(1), final(1)], [intercept(2), final(2)], 'color', 'red'); 

	% Draw tangent and normal
	syms a b c;
	if p.show_tan
		tang = subs(p.TanSlope, {a, b}, {intercept(1), t}) * (p.x-intercept(1)) + intercept(2);
		[~, i] = min(abs(tang-intercept(2)));
		plot(p.x(i-p.ti:i+p.ti), tang(i-p.ti:i+p.ti), '--', 'color', 'black')
	end
	if p.show_norm
		norm = -1/subs(p.TanSlope, {a, b}, {intercept(1), t}) * (p.x-intercept(1)) + intercept(2);
		[~, i] = min(abs(norm-intercept(2)));
		plot(p.x(i-p.ni:i+p.ni), norm(i-p.ni:i+p.ni), '--', 'color', 'black')
	end
end
