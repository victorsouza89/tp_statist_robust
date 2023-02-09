function z = generateTDist(m,v,sigma)
    tau = random('Gamma',v/2,2/v,1,m).^-1;
    n = ((sqrt(sigma))*(random('Normal',0,1,m,1)))';
    z = sqrt(tau).*n;
end