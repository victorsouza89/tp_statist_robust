function z = createTDistribution(v, sigma, N)
    m = size(sigma, 1);
    z = zeros(N, m);
    for i = 1:N
        tau = sqrt(gamrnd(v/2,2/v,1,m).^(-1));
        n = (sqrtm(sigma/2)*(randn(m,1)+1i*randn(m,1)))';
        z(i,:) = sqrt(tau).*n;
    end
end