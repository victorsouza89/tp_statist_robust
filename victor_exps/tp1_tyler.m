clear
close all
clc

m = 2;
r = 0.1;
theta = 2*pi/m;
sigma = get_sigma(m, r, theta);

%% tries to estimate (SCM) sigma for a lot of Ns
v = 10;
Ns = m:m+40;
MC = 100;

err_N = zeros(size(Ns));
err_MC = zeros(length(Ns), MC, m, m);
% monte carlo iteration
for mc = 1:MC
    for k = 1:length(Ns)
        N = Ns(k);
        [z] = generate_tv(v, sigma, N);
        % Tyler estimation starts
        sigma_Tyl_est = estimation_tyler(z, 10^(-3));
        % sigma_Tyl_est is the estimated sigma
        err_MC(k, mc, :, :) = (sigma - sigma_Tyl_est)*(sigma - sigma_Tyl_est)';
    end
end
for k = 1:length(Ns)
    err_MC_i = reshape(sum(err_MC(k,:,:,:)),m,m) * 1/MC;
    err_N(k) = norm(err_MC_i, 'fro');
end

figure
plot(Ns, err_N)
grid on


%% tries to estimate (SCM) sigma for a lot of vs
vs = 0.1:0.1:20;
N = 10;
MC = 100;

err_v = zeros(size(vs));
err_MC = zeros(length(vs), MC, m, m);
% monte carlo iteration
for mc = 1:MC
    for i = 1:length(vs)
        v = vs(i);
        [z] = generate_tv(v, sigma, N);
        % Tyler estimation starts
        sigma_Tyl_est = estimation_tyler(z, 10^(-3));
        % sigma_Tyl_est is the estimated sigma
        err_MC(i, mc, :, :) = (sigma - sigma_Tyl_est)*(sigma - sigma_Tyl_est)';
    end
end
for i = 1:length(vs)
    err_MC_i = reshape(sum(err_MC(i,:,:,:)),m,m) * 1/MC;
    err_v(i) = norm(err_MC_i, 'fro');
end

figure
plot(vs, err_v)
grid on