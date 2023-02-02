function [sigma] = get_sigma(m, r, theta)
    first_vector = ones(1, m);
    for k=2:m
        first_vector(k) = r*exp(1i*theta*(k-1));
    end
    sigma = toeplitz(first_vector);
end