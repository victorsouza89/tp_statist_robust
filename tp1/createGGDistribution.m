function z = createGGDistribution(s, sigma, N)
    m = size(sigma, 1);
    z = zeros(N, m);
    for i = 1:N
        b = (m*gamma(m/s)/gamma((m+1)/s))^s;
        Q = gamrnd(m/s,b,1,m).^(1/s);
        w = ((randn(m,1)+1i*randn(m,1)))';
        u = w/norm(w);
        z(i,:) = sqrt(Q).*sqrtm(sigma)*u';
    end
end