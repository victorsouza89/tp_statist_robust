function e_rms = estimatorSCM(n,m,v,sigma)
    z = zeros(n,m);
    for k = 1:n
        z(k,:) = generateTDist(m,v,sigma);
    end
    sigmaSCM = zeros(m,m);
    for k = 1:n
        sigmaSCM = sigmaSCM + z(k,:)*z(k,:)';
    end
    sigmaSCM = sigmaSCM/n;
    %sigmaSCM = sigmaSCM/trace(sigmaSCM);
    sigmav = reshape(sigma,1,m*m);
    sigmaSCMv = reshape(sigmaSCM,1,m*m);
    e_rms = norm(mean((sigmav-sigmaSCMv)*(sigmav-sigmaSCMv)'),'fro');
end