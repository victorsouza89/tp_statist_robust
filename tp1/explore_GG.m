%% variation with n
clear variables; clc;
close all; 

% model parameters
m = 2;
r = 0.1;
theta = 2*pi/m;
sigma = get_sigma(m, r, theta);
get_err = @(sigmaEST) norm(reshape(sigma-sigmaEST,m*m,1)'*reshape(sigma-sigmaEST,m*m,1),'fro');

% huber estimator parameter 
q = 0.5;

% simulation parameters
N = 40;
s = 10;
MC = 100;

% init arrays
ERMS = zeros(1,length(N));
ERMS_tyler = zeros(1, length(N));
ERMS_huber = zeros(1, length(N));

% simulation
for n = 1:N
    eMC = zeros(1,MC);
    eMCTY = zeros(1,MC);
    eMCHU = zeros(1,MC);
    % Monte Carlo
    for k = 1:MC
        % create data
        z = createGGDistribution(s,sigma,n);
        % estimates
        sigmaSCM = z'*z/n;
        sigmaCSCM = m*sigmaSCM/trace(sigmaSCM);
        sigmaTYLER = calculateTylerEstimator(m,n,z);
        sigmaHUBER = calculateHuberEstimator(m,n,q,z);
        % saves estimation
        eMC(k) = get_err(sigmaCSCM);
        eMCTY(k) = get_err(sigmaTYLER);
        eMCHU(k) = get_err(sigmaHUBER);
    end
    ERMS(n) = mean(eMC);
    ERMS_tyler(n) = mean(eMCTY);
    ERMS_huber(n) = mean(eMCHU);
end


%% Print results in function of n
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

%% variation with s
%clear variables; clc;
%close all; 

% model parameters
m = 2;
r = 0.1;
theta = 2*pi/m;
sigma = get_sigma(m, r, theta);
get_err = @(sigmaEST) norm(reshape(sigma-sigmaEST,m*m,1)'*reshape(sigma-sigmaEST,m*m,1),'fro');

% huber estimator parameter 
q = 0.5;

% simulation parameters
S = 0.1:0.1:10;
n = 10;
MC = 100;

% init arrays
ERMS = zeros(1,length(S));
ERMS_tyler = zeros(1, length(S));
ERMS_huber = zeros(1, length(S));

% simulation
for s_aux = 1:length(S)
    s = S(s_aux);
    eMC = zeros(1,MC);
    eMCTY = zeros(1,MC);
    eMCHU = zeros(1,MC);
    % Monte Carlo
    for k = 1:MC
        % create data
        z = createGGDistribution(s,sigma,n);  
        % estimates
        sigmaSCM = z'*z/n;
        sigmaCSCM = m*sigmaSCM/trace(sigmaSCM);
        sigmaTYLER = calculateTylerEstimator(m,n,z);
        sigmaHUBER = calculateHuberEstimator(m,n,q,z);
        % saves estimation
        eMC(k) = get_err(sigmaCSCM);
        eMCTY(k) = get_err(sigmaTYLER);
        eMCHU(k) = get_err(sigmaHUBER);
    end
    ERMS(s_aux) = mean(eMC);
    ERMS_tyler(s_aux) = mean(eMCTY);
    ERMS_huber(s_aux) = mean(eMCHU);
end
%% Print results in function of s
figure 
plot(S,ERMS)
hold on 
plot(S,ERMS_tyler)
hold on 
plot(S,ERMS_huber)
grid on
xlabel('s (degree of freedom)')
ylabel('ERMS value')
legend('SCM','Tyler','Huber')