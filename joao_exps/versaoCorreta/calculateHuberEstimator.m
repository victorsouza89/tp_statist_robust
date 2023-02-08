function sigma_huber = calculateHuberEstimator(m,n,c,z)
    q = chi2cdf(2*c^2,2*m);
    b = chi2cdf(2*c^2,2*m+2) + c^2*(1-q)/m;
    tol = 10^-3;
    sigma_huberk = eye(m);
    sigma_huber = zeros(m);
    while norm(sigma_huber-sigma_huberk,'fro') > norm(sigma_huber)*tol
        S = zeros(m);
        for i = 1:n
            zi = z(i,:);
            tk = zi*inv(sigma_huberk)*zi';
            if tk <= c^2
                phi_huber = 1/b;
            else
                phi_huber = c^2/(tk*b);
            end
            S = S + phi_huber.*z'*z;
        end
        S = S/n;
        sigma_huberk = m*S/trace(S);
        sigma_huber = sigma_huberk;
    end
end