function sigma = calculateHuberEstimator(m,n,q,z)
    c = sqrt(chi2inv(q,2*m)/2);
    b = chi2cdf(2*c^2,2*m+2) + c^2*(1-q)/m;
    tol = 10^-2;
    sigma_k = eye(m);
    sigma = zeros(m);
    while norm(sigma-sigma_k,'fro') > tol * norm(sigma,'fro')
        S = zeros(m);
        for i = 1:n
            zi = z(i,:);
            tk = zi/sigma_k*zi';
            % phi calc start
            if tk <= c^2
                phi = 1/b;
            else
                phi = c^2/(tk*b);
            end
            % phi calc end
            S = S + phi.*zi'*zi;
        end
        S = S/n;
        sigma_k = m*S/trace(S);
        sigma = sigma_k;
    end
end