function [z,sigma] = createTDistribution(n,m,v,r)
    z = zeros(n,m);
    theta = 2*pi/m;
    vec_toeplitz = ones(1,m);
    for k = 2:m
        vec_toeplitz(k) = r*exp(1i*theta*(k-1));
    end
    sigma = toeplitz(vec_toeplitz);
    for k = 1:n    
        tau = sqrt(gamrnd(v/2,2/v,1,m).^(-1));
        n = (sqrtm(sigma/2)*(randn(m,1)+1i*randn(m,1)))';
        z(k,:) = sqrt(tau).*n;
    end

    
end