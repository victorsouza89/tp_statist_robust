function sigma = calculateTylerEstimator(m,n,z)
    tol = 10^-2;
    sigma_k = eye(m);
    sigma = zeros(m);
    while norm(sigma-sigma_k,'fro') > tol * norm(sigma,'fro')
        S = zeros(m);
        for i = 1:n
            zi = z(i,:);
            tk = zi/sigma_k*zi';
            % phi calc start
            phi = m/tk;
            % phi calc end
            S = S + phi.*zi'*zi;
        end
        S = S/n;
        sigma_k = m*S/trace(S);
        sigma = sigma_k;
    end
    sigma = 0.5*sigma;
end