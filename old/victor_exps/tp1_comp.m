clear
close all
clc

m = 5;
r = 0.1;
theta = 2*pi/m;
sigma = get_sigma(m, r, theta);

%% tries to estimate (SCM) sigma for a lot of Ns
v = 10;
Ns = m:m+40;
MC = 100;

err_Tyl_N = zeros(size(Ns));
err_SCM_N = zeros(size(Ns));
err_Tyl_MC = zeros(length(Ns), MC, m, m);
err_SCM_MC = zeros(length(Ns), MC, m, m);
% monte carlo iteration
for mc = 1:MC
    for k = 1:length(Ns)
        N = Ns(k);
        [z] = generate_tv(v, sigma, N);
        % Tyler estimation starts
        sigma_Tyl_est = estimation_tyler(z, 10^(-3));
        sigma_SCM_est = estimation_SCM(z);
        % sigma_Tyl_est is the estimated sigma
        err_Tyl_MC(k, mc, :, :) = (sigma - sigma_Tyl_est)*(sigma - sigma_Tyl_est)';
        err_SCM_MC(k, mc, :, :) = (sigma - sigma_SCM_est)*(sigma - sigma_SCM_est)';
    end
end
for k = 1:length(Ns)
    err_Tyl_N(k) = norm(sum(err_Tyl_MC(k,:,:,:)) * 1/MC, 'fro');
    err_SCM_N(k) = norm(sum(err_SCM_MC(k,:,:,:)) * 1/MC, 'fro');
end

figure
plot(Ns, err_Tyl_N)
hold on
plot(Ns, err_SCM_N)
grid on
xlabel('N') 
ylabel('RMSE') 



%% tries to estimate (SCM) sigma for a lot of vs
vs = 0.1:0.1:20;
N = 10;
MC = 100;

err_Tyl_v = zeros(size(vs));
err_SCM_v = zeros(size(vs));
err_Tyl_MC = zeros(length(vs), MC, m, m);
err_SCM_MC = zeros(length(vs), MC, m, m);
% monte carlo iteration
for mc = 1:MC
    for i = 1:length(vs)
        v = vs(i);
        [z] = generate_tv(v, sigma, N);
        % Tyler estimation starts
        sigma_Tyl_est = estimation_tyler(z, 10^(-5));
        sigma_SCM_est = estimation_SCM(z);
        % sigma_Tyl_est is the estimated sigma
        err_Tyl_MC(k, mc, :, :) = (sigma - sigma_Tyl_est)*(sigma - sigma_Tyl_est)';
        err_SCM_MC(k, mc, :, :) = (sigma - sigma_SCM_est)*(sigma - sigma_SCM_est)';
    end
end
for k = 1:length(Ns)
    err_Tyl_v(k) = norm(sum(err_Tyl_MC(k,:,:,:)) * 1/MC, 'fro');
    err_SCM_v(k) = norm(sum(err_SCM_MC(k,:,:,:)) * 1/MC, 'fro');
end

figure
plot(vs, err_Tyl_v)
hold on
plot(vs, err_SCM_v)
grid on
xlabel('v') 
ylabel('RMSE') 
