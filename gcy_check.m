clear;
clc;
format long
tic;

%% paramter

L = 1220;
lambda = 3;
% beta_all = [0.1:0.1:10-0.1,10:1:100-1,100:10:1000-10];
beta_all = 1+0.01:0.02:10;
nbeta = length(beta_all);
alpha = 377/610;

H1 = zeros(L);
for i = 1:L-1
    H1(i,i+1) = 1;
    H1(i+1,i) = 1;
end
H1(L,1) = 1;
H1(1,L) = 1;

phi4 = zeros(nbeta,1);

for n = 1:nbeta
    beta = beta_all(n);
    H2 = zeros(L,1);
    for i = 1:L/2
        H2(2*i) = -tanh(beta*(cos(2*pi*alpha*i)-cos(pi*alpha)))/tanh(beta);
    end

    H = H1 + lambda*diag(H2);
    [V, D] = eig(H);
%     e = diag(D);
%     [V, D] = eigs(sparse(H),2,'smallestreal','Tolerance',1e-4);
    phi4(n) = sum(abs(V(:,1)).^4);
end

figure
plot(beta_all,phi4);
set(gca,'XScale','log')

toc;