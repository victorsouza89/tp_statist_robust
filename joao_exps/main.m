clear variables; close all;
%% variation with n
m = 2;
N = 40;
v = 10;
r = 0.1;
%z = createTDistribution(N,m,v,r); generating N variables z of length m

MC = 100;
ERMS = zeros(1,length(N));

for n = 1:N
    eMC = zeros(1,MC);
    for k = 1:MC
        sigmaSCM = zeros(m,m);
        sigmaCSCM = zeros(m,m);
        [z,sigma] = createTDistribution(n,m,v,r);  
        sigmaSCM = z'*z/n;
        sigmaCSCM = m*sigmaSCM/trace(sigmaSCM);
        eMC(k) = norm(reshape(sigma-sigmaCSCM,m*m,1)'*reshape(sigma-sigmaCSCM,m*m,1),'fro');    
    end
    ERMS(n) = mean(eMC);
end
%%
figure 
plot(1:N,ERMS)
grid on
xlabel('N (number of z variables)')
ylabel('ERMS value')

%% variation with v 
clear variables; close all;
V = 0.1:0.1:10;
n = 10;
m = 2;
r = 0.1;
MC = 100;

for v_aux = 1:length(V)
    v = V(v_aux);
    eMC = zeros(1,MC);
    for k = 1:MC
        sigmaSCM = zeros(m,m);
        sigmaCSCM = zeros(m,m);
        [z,sigma] = createTDistribution(n,m,v,r);  
        sigmaSCM = z'*z/n;
        sigmaCSCM = m*sigmaSCM/trace(sigmaSCM);
        eMC(k) = norm(reshape(sigma-sigmaCSCM,m*m,1)'*reshape(sigma-sigmaCSCM,m*m,1),'fro');    
    end
    ERMS(v_aux) = mean(eMC);
end
%%
figure 
plot(V,ERMS)
grid on
xlabel('v (degree of freedom)')
ylabel('ERMS value')