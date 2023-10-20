clear;
clc;
format long
tic;

%% paramter
myseed = 1;
rng(myseed)

L = 10;
num_T = 2000;
T = 0:2*num_T;
nT = 2*num_T+1;
len = 2^L;
g0 = 0.4;
J0 = -1;
t = 2;
% g = g0*(1+0.01*randn(L,1));
% J = J0*(1+1*randn(L,1));
% g = g0*(2*round(rand(L,1))-1);
% J = J0*(2*round(rand(L,1))-1);
g = g0*ones(L,1);
J = J0*ones(L,1);

sigmaz = [1;-1];
sigmax = [0 1;1 0];
I2 = eye(2);
II2 = ones(2);

%% construction of Hamiltonian and observable
Hz = zeros(len,1);
for i = 1:L-1
    Hz_temp = ones(2^(i-1),1);
    Hz_temp = kron(Hz_temp, sigmaz);
    Hz_temp = kron(Hz_temp, sigmaz);
    Hz_temp = kron(Hz_temp, ones(2^(L-i-1),1));
    Hz = Hz + J(i)*Hz_temp;
end
% PBC
Hz_temp = sigmaz;
Hz_temp = kron(Hz_temp, ones(2^(L-2),1));
Hz_temp = kron(Hz_temp, sigmaz);
Hz = Hz + J(L)*Hz_temp;

Hx = zeros(len);
matrix_sz = zeros(len,L);
for i = 1:L
    Hx_temp = eye(2^(i-1));
    Hx_temp = kron(Hx_temp,sigmax);
    Hx_temp = kron(Hx_temp,eye(2^(L-i)));
    Hx = Hx + g(i)*Hx_temp;

    sz_temp = ones(2^(i-1),1);
    sz_temp = kron(sz_temp,sigmaz);
    sz_temp = kron(sz_temp,ones(2^(L-i),1));
    matrix_sz(:,i) = sz_temp;
end

%% initial state

% phi0 = rand(len,1);
% phi0 = phi0./sqrt(sum(phi0.^2));

[~,index0] = min(Hz);
phi0 = zeros(len,1);
pos = index0(1);

sz_GS0 = zeros(L,1);
temp = pos-1;
for i = 1:L
    sz_GS0(i) = mod(temp,2);
    temp = floor(temp/2);
end
sz_GS = 2*sz_GS0-1;
phi0(pos) = 1;

% % flip a single point
% sz_GS_new = sz_GS;
% sz_GS_new(floor(L/2)) = - sz_GS(floor(L/2));
% sz_GS0_new = (sz_GS_new+1)/2;
% temp = 0;
% for i = 1:L
%     temp = temp*2;
%     temp = temp + sz_GS0_new(L-i+1);    
% end
% temp = temp+1;
% phi0(temp) = 1;
% phi0 = phi0./sqrt(sum(phi0.^2));

phi = phi0;
phit_store = zeros(len,nT);
phit_store(:,1) = phi0;

%% time evolution

H1 = diag(Hz) + Hx;
H2 = diag(Hz) - Hx;
[V1,D1] = eig(H1);
e1 = diag(D1);
[V2,D2] = eig(H2);
e2 = diag(D2);

trans1 = exp(-1i*t*e1);
trans2 = exp(-1i*t*e2);

now = 2;
for i = 1:num_T
    % H1
    temp = V1'*phi;
    temp = temp.*trans1;
    phi = V1*temp;
    phit_store(:,now) = phi;
    now = now + 1;

    % H2
    temp = V2'*phi;
    temp = temp.*trans2;
    phi = V2*temp;
    phit_store(:,now) = phi;
    now = now + 1;
end

%% calculate observable

sz = zeros(L,nT);
for i = 1:L
    sz(i,:) = sum(conj(phit_store).*(matrix_sz(:,i).*phit_store));
end

sz_mean = sum(sz_GS.*sz)/L;

figure;
set(gcf, 'position', [250 70 1400 900]);
subplot(2,1,1)
plot(T,sz_mean);
subplot(2,1,2)
imagesc(sz);

toc;
