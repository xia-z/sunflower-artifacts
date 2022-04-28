%Simulator: simution of using fft for decoding ook signal under different frequencies

close all;
clearvars;
clear;

fs = 5*1e4;  % samplong rate

f1 = 1*1e4; % 10k
f3 = 500; % frequecy 1 for fsk
f4 = 1000;% frequecy 2 for fsk

t0 = 0:1/fs:0.0005-1/fs;
v1 = (square(2*pi*t0*f1)+1)/2; % 1  ook   needs to be 0.0005 s
v2 = 0.8*((square(2*pi*t0*f1)+1)/2); % 0  ook  0.0005 s


% v3 = [v1,v1,v1,v1,v1,v1,v1,v1,v1,v1,v2,v2,v2,v2,v2,v2,v2,v2,v2,v2]; %fsk 500  0
% v4 = [v1,v2,v1,v2,v1,v2,v1,v2,v1,v2,v1,v2,v1,v2,v1,v2,v1,v2,v1,v2,]; %fsk 5k   1


v3 = [v1,v1,v2,v2]; %fsk 500  0
v4 = [v1,v2,v1,v2]; %fsk 1k   1

% dv1 = [diff(v1),1];
% dv2 = [diff(v2),1];
% dv3 = [diff(v3),1];

v = [v3,v4,v4,v3,v3,v4,v3]; % 0110010
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



n = 100;  % number of FFT.  determines the  f-resolution = fs/n  we want it to be 500 at least

a = buffer(x,n,n-1); % overlap . padding in the beginning but not the end

i = 1;
y_want = [];
coe = [];

while i < length(a(1,:))+1

%     R = corrcoef(a(:,i), dv3');
    Y = fft(a(:,i), n);
    P = abs(Y/n).^2;
    P1= P(1:n/2+1);
    y_want = [y_want; P1(2), P1(3)];
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
hold off
% plot(1:length(y_want), y_want(:,3), 'r', 'LineWidth',3);
% hold off;

% legend('200', '400', '500');
legend('0.5k', '1k');
xlabel("Sliding Window ");
ylabel("Normalized Power");

