clear
close all
clc

m = 2;
r = 0.1;
theta = 2*pi/m;
sigma = get_sigma(m, r, theta);

%% tries to estimate (SCM) sigma for a lot of Ns
v = 10;
N = 40;
MC = 100;

err_N = zeros(1,N);
err_MC = zeros(N, MC, m, m);
% monte carlo iteration
for mc = 1:MC
    sigma_SCM = zeros(N, m, m);
    % SCM estimation starts
    for n = 1:N
        [z] = generate_tv(v, sigma, n);
        sigma_SCM(n, :, :) = (z'*z);
        sigma_SCM_est = (1/n) * reshape(sum(sigma_SCM),m,m);
        sigma_SCM_est = sigma_SCM_est * m/trace(sigma_SCM_est);
        % sigma_SCM_est is the estimated sigma (for this n)
        err_MC(n, mc, :, :) = (sigma - sigma_SCM_est)*(sigma - sigma_SCM_est)';
    end
end
for n = 1:N
    err_MC_i = reshape(sum(err_MC(n,:,:,:)),m,m) * 1/MC;
    err_N(n) = norm(err_MC_i, 'fro');
end

figure
plot(1:N, err_N)
grid on


%% tries to estimate (SCM) sigma for a lot of vs
vs = 0.1:0.1:10;
N = 10;
MC = 1000;

err_v = zeros(size(vs));
err_MC = zeros(length(vs), MC, m, m);
% monte carlo iteration
for mc = 1:MC
    for k = 1:length(vs)
        v = vs(k);
        [z] = generate_tv(v, sigma, N);
        % SCM estimation starts
        sigma_SCM_est = estimation_SCM(z);
        % sigma_SCM_est is the estimated sigma
        err_MC(k, mc, :, :) = (sigma - sigma_SCM_est)*(sigma - sigma_SCM_est)';
    end
end
for k = 1:length(vs)
    err_MC_i = reshape(sum(err_MC(k,:,:,:)),m,m) * 1/MC;
    err_v(k) = norm(err_MC_i, 'fro');
end

figure
plot(vs, err_v)
grid on