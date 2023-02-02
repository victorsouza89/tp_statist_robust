function sigma = generateSigmaTDist(r,m)
    theta = 2*pi/m;
    vec_toeplitz = ones(1,m);
    for k = 2:m
        vec_toeplitz(k) = r*exp(1i*theta*(k-1));
    end
    sigma = toeplitz(vec_toeplitz);
end
