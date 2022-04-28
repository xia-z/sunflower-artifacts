%Simulator: simution of using fft for decoding fsk signal(high frequecies) under different frequencies(lower frequecies)

close all;
clearvars;
clear;

fs = 10*1e4; % samplong rate

f1 = 1*1e4; %10 k     frequecy 1 for higher frequency pair
f2 = 2*1e4; %20 k    frequecy 2 for higher frequency pair

f3 = 500; %  frequecy 1 for lower frequency pair
f4 = 1000;%  frequecy 1 for lower frequency pair

t0 = 0:1/fs:0.0005-1/fs;
v1 = (square(2*pi*t0*f1)+1)/2; % 1  ook   needs to be 0.0005 s
v2 = 0.5*((square(2*pi*t0*f1)+1)/2); % 0  ook  0.0005 s

v21 = (square(2*pi*t0*f2)+1)/2; % 1  ook   needs to be 0.0005 s
v22 = 0.5*((square(2*pi*t0*f2)+1)/2); % 0  ook  0.0005 s

% v3 = [v1,v1,v1,v1,v1,v1,v1,v1,v1,v1,v2,v2,v2,v2,v2,v2,v2,v2,v2,v2]; %fsk 500  0
% v4 = [v1,v2,v1,v2,v1,v2,v1,v2,v1,v2,v1,v2,v1,v2,v1,v2,v1,v2,v1,v2,]; %fsk 5k   1

%fs 

v3 = [v1,v1,v2,v2]; %fsk 500  0
v4 = [v1,v2,v1,v2]; %fsk 1k   1

v5 = [v21,v21,v22,v22]; %fsk 500  0
v6 = [v21,v22,v21,v22]; %fsk 1k   1

% dv1 = [diff(v1),1];
% dv2 = [diff(v2),1];
% dv3 = [diff(v3),1];

v = [v3,v6,v4,v5,v3,v6,v4,v5]; % fsk1: 01100110  %fsk2: 0 1 0 1 0 0101
% dv = diff(v);
% dv = [dv,1];
% x = v;
t = 0:1/fs:0.0005*length(v)/length(v1)-1/fs;
x = v;
% t = t(2:end);

figure(1)
plot(t, x);
xlabel("Time [s]")
ylabel("Signal [V]")



n = 200;  % number of FFT.  determines the  f-resolution = fs/n  we want it to be 500 at least

a = buffer(x,n,n-1); % overlap . padding in the beginning but not the end

i = 1;
y_want = [];
coe = [];

while i < length(a(1,:))+1

%     R = corrcoef(a(:,i), dv3');
    Y = fft(a(:,i), n);
    P = abs(Y/n).^2;
    P1= P(1:n/2+1);
    y_want = [y_want; P1(2), P1(3), P1(21), P1(41)];
    i = i +1;
end

f_plot = fs*(0:(n/2))/n;




figure(2)
plot(f_plot, P1);
xlabel("Frequecies ")
ylabel("Normalized Power")

plot(1:length(y_want), y_want(:,1), 'b', 'LineWidth',3);
hold on
plot(1:length(y_want), y_want(:,2), 'LineWidth',3);
hold on
plot(1:length(y_want), y_want(:,3), 'r', 'LineWidth',3);
hold on
plot(1:length(y_want), y_want(:,4), 'y','LineWidth',3);
hold off;

% legend('200', '400', '500');
legend('0.5k', '1k','10k','20k');
xlabel("Sliding Window ");
ylabel("Normalized Power");

