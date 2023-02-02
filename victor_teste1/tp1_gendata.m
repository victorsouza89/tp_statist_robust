clear
close all
clc

v = 10;

%% generates pdf t_v
step = 0.1;
x =-40:step:40;
y = cdf('T',x,v);
dev_y = zeros(1,length(y));
for k = 2 : length(y)
    dev_y(k) = (y(k) - y(k-1))/step;
end
figure
plot(x,dev_y)

%% generates data t_v
m = 10;
r = 0.001;
theta = 2*pi/m;
sigma = get_sigma(m, r, theta);
[z] = generate_tv(v, sigma, 1);

figure
histogram(real(z), [-40:1:40])
z_ = random('T', v, 1, m);
figure
histogram(real(z_), [-40:1:40])

