function sigma_tyler = calculateTylerEstimator(m,n,z)
    tol = 10^-2;
    sigma_tylerk = eye(m);
    sigma_tyler = zeros(m);
    while norm(sigma_tyler-sigma_tylerk,'fro') > tol
        S = zeros(m);
        for i = 1:n
            zi = z(i,:);
            tk = zi*inv(sigma_tylerk)*zi';
            phi_tyler = m/tk;
            S = S + phi_tyler.*zi'*zi;
        end
        S = S/n;
        sigma_tylerk = m*S/trace(S);
        sigma_tyler = sigma_tylerk;
    end
    sigma_tyler = 0.5*sigma_tyler;
end