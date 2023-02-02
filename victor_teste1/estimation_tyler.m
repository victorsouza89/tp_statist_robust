function [sigma_Tyl_est] = estimation_tyler(z, tol)
    N = size(z, 1);
    m = size(z, 2);
    sigma_Tyl = zeros(N, m, m);
    sigma_Tyl_est = eye(m);
    sigma_Tyl_est_old = zeros(m, m);
    while norm(sigma_Tyl_est - sigma_Tyl_est_old, 'fro') > norm(sigma_Tyl_est, 'fro') * tol 
        for n = 1:N
            t_k = z(n,:) / sigma_Tyl_est * z(n,:)';
            phi_k = m / t_k;
            sigma_Tyl(n, :, :) = phi_k * (z'*z);
        end
        sigma_Tyl_est = (1/N) * reshape(sum(sigma_Tyl),m,m);
        sigma_Tyl_est = sigma_Tyl_est * m/trace(sigma_Tyl_est);
        sigma_Tyl_est_old = sigma_Tyl_est;
    end
end