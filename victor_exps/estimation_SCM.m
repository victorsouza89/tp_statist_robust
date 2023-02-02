function [sigma_SCM_est] = estimation_SCM(z)
    N = size(z, 1);
    m = size(z, 2);
    sigma_SCM = zeros(N, m, m);
    for n = 1:N
            sigma_SCM(n, :, :) = (z'*z);
    end
    sigma_SCM_est = (1/N) * reshape(sum(sigma_SCM),m,m);
    sigma_SCM_est = sigma_SCM_est *m/trace(sigma_SCM_est);
end