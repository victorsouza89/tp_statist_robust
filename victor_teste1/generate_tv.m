function [z, sigma] = generate_tv(v, sigma, N)
    m = size(sigma, 1);
    
    z = zeros(N, m);
    for i = 1:N
        tau_1 = random("Gamma", v/2, 2/v, 1, m);
        tau = 1./tau_1;
        n = sqrt(sigma/2) * (randn(m,1) + 1i*randn(m,1));
        n = n.';
        z(i,:) = sqrt(tau) .* n;
    end
end