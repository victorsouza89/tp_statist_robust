clear variables
close all 
%% Generating z (t-distributed data)
% m = 1000;
% v = 1;
% tau = random('Gamma',v/2,2/v,1,m).^(-1);
% r = 0.5;
% theta = 2*pi/m;
% vec_toeplitz = ones(1,m);
% for k = 2:m
%     vec_toeplitz(k) = r*exp(1i*theta*(k-1));
% end
% sigma = toeplitz(vec_toeplitz);
% n = (sqrt(sigma/2))*(randn(m,1)+1i*randn(m,1));
% n = n.';
% z = sqrt(tau).*n;
m = 1000;
% Tau
figure
sgtitle("Tau")
subplot(2,1,1)
v = 10;
tau = random('Gamma',v/2,2/v,1,m).^-1;
grid on
plot(tau)
subplot(2,1,2)
histogram(tau)
grid on 

% n zero mean complex Gaussian random 
r = 0.0001;
theta = 2*pi/m;
vec_toeplitz = ones(1,m);
for k = 2:m
    vec_toeplitz(k) = r*exp(1i*theta*(k-1));
end
sigma = toeplitz(vec_toeplitz);
n = ((sqrt(sigma))*(random('Normal',0,1,m,1)))';
figure
sgtitle("CN (real part)")
subplot(2,1,1)
grid on
plot(real(n))
subplot(2,1,2)
histogram(real(n))
grid on 

% z (t-distributed function)
%z = sqrt(tau).*n;
z = generateTDist(m,v,r);
figure
sgtitle("t-distributed data")
subplot(2,1,1)
grid on
plot(real(z))
subplot(2,1,2)
histogram(real(z))
grid on 

%% Generating a set of n t-distributed data
n = 20;
m = 10;
z = zeros(n,m);
v = 10;
r = 0.5;
% figure
for k = 1:n
    z(k,:) = generateTDist(m,v,r);
%     plot(real(z(k,:)))
%     hold on
end

%% SCM Estimator
sigmaSCM = zeros(m,m);
% figure
for k = 1:n
    sigmaSCM = sigmaSCM + z(k,:)*z(k,:)';  
%     plot(real(z(k,:)))
%     hold on
end
sigmaSCM = sigmaSCM/n;
%sigmaSCM = sigmaSCM/trace(sigmaSCM);
sigma = generateSigmaTDist(r,m);
sigmav = reshape(sigma,1,m*m);
sigmaSCMv = reshape(sigmaSCM,1,m*m);
E_SCM = norm(mean((sigmav-sigmaSCMv)*(sigmav-sigmaSCMv)'),'fro');

%%
m = 10;
sigma = generateSigmaTDist(r,m);
n = m:1:1000;
e_scm = zeros(1,length(n));
for k =1:length(n)    
    e_scm(k) = estimatorSCM(n(k),m,v,sigma);
end
figure
plot(n,e_scm)
xlabel('n')
ylabel('ERMS')