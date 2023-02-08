clear variables; close all;
%% variation with n
m = 2;
N = 40;
v = 10;
r = 0.1;
%z = createTDistribution(N,m,v,r); generating N variables z of length m

MC = 100;
ERMS = zeros(1,length(N));
ERMS_tyler = zeros(1, length(N));
ERMS_huber = zeros(1, length(N));
for n = 1:N
    eMC = zeros(1,MC);
    eMCTY = zeros(1,MC);
    eMCHU = zeros(1,MC);
    for k = 1:MC
        sigmaSCM = zeros(m,m);
        sigmaCSCM = zeros(m,m);
        sigmaTYLER = zeros(m,m);
        sigmaHUBER = zeros(m);
        [z,sigma] = createTDistribution(n,m,v,r);  
        sigmaSCM = z'*z/n;
        sigmaCSCM = m*sigmaSCM/trace(sigmaSCM);
        sigmaTYLER = calculateTylerEstimator(m,n,z);
        sigmaHUBER = calculateHuberEstimator(m,n,1,z);
        eMC(k) = norm(reshape(sigma-sigmaCSCM,m*m,1)'*reshape(sigma-sigmaCSCM,m*m,1),'fro');    
        eMCTY(k) = norm(reshape(sigma-sigmaTYLER,m*m,1)'*reshape(sigma-sigmaTYLER,m*m,1),'fro');
        eMCHU(k) = norm(reshape(sigma-sigmaHUBER,m*m,1)'*reshape(sigma-sigmaHUBER,m*m,1),'fro');
    end
    ERMS(n) = mean(eMC);
    ERMS_tyler(n) = mean(eMCTY);
    ERMS_huber(n) = mean(eMCHU);
end



%% Results in function of n
figure 
plot(1:N,ERMS)
hold on
plot(1:N,ERMS_tyler)
hold on 
plot(1:N,ERMS_huber)
grid on
xlabel('N (number of z variables)')
ylabel('ERMS value')
legend('SCM','Tyler','Huber')

%% variation with v 
clear variables; close all;
V = 0.1:0.1:10;
n = 10;
m = 2;
r = 0.1;
MC = 100;
ERMS = zeros(1,length(V));
ERMS_tyler = zeros(1, length(V));
ERMS_huber = zeros(1, length(V));

for v_aux = 1:length(V)
    v = V(v_aux);
    eMC = zeros(1,MC);
    eMCTY = zeros(1,MC);
    eMCHU = zeros(1,MC);
    for k = 1:MC
        sigmaSCM = zeros(m,m);
        sigmaCSCM = zeros(m,m);
        sigmaTYLER = zeros(m);
        sigmaHUBER = zeros(m);
        [z,sigma] = createTDistribution(n,m,v,r);  
        sigmaSCM = z'*z/n;
        sigmaCSCM = m*sigmaSCM/trace(sigmaSCM);
        sigmaTYLER = calculateTylerEstimator(m,n,z);
        sigmaHUBER = calculateHuberEstimator(m,n,1,z);
        eMC(k) = norm(reshape(sigma-sigmaCSCM,m*m,1)'*reshape(sigma-sigmaCSCM,m*m,1),'fro');  
        eMCTY(k) = norm(reshape(sigma-sigmaTYLER,m*m,1)'*reshape(sigma-sigmaTYLER,m*m,1),'fro');
        eMCHU(k) = norm(reshape(sigma-sigmaHUBER,m*m,1)'*reshape(sigma-sigmaHUBER,m*m,1),'fro');
    end
    ERMS(v_aux) = mean(eMC);
    ERMS_tyler(v_aux) = mean(eMCTY);
    ERMS_huber(n) = mean(eMCHU);
end
%%
figure 
plot(V,ERMS)
hold on 
plot(V,ERMS_tyler)
hold on 
plot(V,ERMS_huber)
grid on
xlabel('v (degree of freedom)')
ylabel('ERMS value')
legend('SCM','Tyler','Huber')