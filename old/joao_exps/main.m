clear variables; close all;
%% variation with n
m = 5;
N = 40;
v = 10;
r = 0.1;
q = 0.01;
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
        sigmaHUBER = calculateHuberEstimator(m,n,q,z);
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
legend('SCM','Tyler','Huber q = '+string(q))
title('t-distributed data \nu = ' +string(v))

%% variation with v 
clear variables; close all;
V = 0.1:0.1:10;
n = 4;
m = 5;
r = 0.1;
q = 0.9;
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
        sigmaHUBER = calculateHuberEstimator(m,n,q,z);
        eMC(k) = norm(reshape(sigma-sigmaCSCM,m*m,1)'*reshape(sigma-sigmaCSCM,m*m,1),'fro');  
        eMCTY(k) = norm(reshape(sigma-sigmaTYLER,m*m,1)'*reshape(sigma-sigmaTYLER,m*m,1),'fro');
        eMCHU(k) = norm(reshape(sigma-sigmaHUBER,m*m,1)'*reshape(sigma-sigmaHUBER,m*m,1),'fro');
    end
    ERMS(v_aux) = mean(eMC);
    ERMS_tyler(v_aux) = mean(eMCTY);
    ERMS_huber(v_aux) = mean(eMCHU);
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
legend('SCM','Tyler','Huber q = '+string(q))
%% GG distributed Data
clear variables; close all;
%%
m = 2 ; 
s = 1 ; 
N = 40;
r = 0.1;
q = 0.5;
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
        sigma = get_sigma(m,r,2*pi/m);
        z = createGGDistribution(s, sigma, n);  
        sigmaSCM = z'*z/n;
        sigmaCSCM = m*sigmaSCM/trace(sigmaSCM);
        sigmaTYLER = calculateTylerEstimator(m,n,z);
        sigmaHUBER = calculateHuberEstimator(m,n,q,z);
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
legend('SCM','Tyler','Huber q = '+string(q))
title('GG-distributed Data s = ' + string(s))

%% variation with s parameter 
clear variables; close all;
S = 0.1:0.1:10;
n = 10;
m = 2;
r = 0.1;
q = 0.5;
MC = 100;
ERMS = zeros(1,length(S));
ERMS_tyler = zeros(1, length(S));
ERMS_huber = zeros(1, length(S));

for s_aux = 1:length(S)
    s = S(s_aux);
    eMC = zeros(1,MC);
    eMCTY = zeros(1,MC);
    eMCHU = zeros(1,MC);
    for k = 1:MC
        sigmaSCM = zeros(m,m);
        sigmaCSCM = zeros(m,m);
        sigmaTYLER = zeros(m);
        sigmaHUBER = zeros(m);
        sigma = get_sigma(m, r, 2*pi/m);
        z = createGGDistribution(s, sigma, n);
        sigmaSCM = z'*z/n;
        sigmaCSCM = m*sigmaSCM/trace(sigmaSCM);
        sigmaTYLER = calculateTylerEstimator(m,n,z);
        sigmaHUBER = calculateHuberEstimator(m,n,q,z);
        eMC(k) = norm(reshape(sigma-sigmaCSCM,m*m,1)'*reshape(sigma-sigmaCSCM,m*m,1),'fro');  
        eMCTY(k) = norm(reshape(sigma-sigmaTYLER,m*m,1)'*reshape(sigma-sigmaTYLER,m*m,1),'fro');
        eMCHU(k) = norm(reshape(sigma-sigmaHUBER,m*m,1)'*reshape(sigma-sigmaHUBER,m*m,1),'fro');
    end
    ERMS(s_aux) = mean(eMC);
    ERMS_tyler(s_aux) = mean(eMCTY);
    ERMS_huber(s_aux) = mean(eMCHU);
end
%%
figure 
plot(S,ERMS)
hold on 
plot(S,ERMS_tyler)
hold on 
plot(S,ERMS_huber)
grid on
xlabel('s (degree of freedom)')
ylabel('ERMS value')
legend('SCM','Tyler','Huber q = '+string(q))
title('GG-distributed data n = '+string(n))