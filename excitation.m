clear;
clc;
format long
tic;

%% paramter
myseed = 1;
rng(myseed)

L = 4;
num_T = 100;
T = 0:2*num_T;
nT = 2*num_T+1;
len = 2^L;
g0 = pi/2+0.1;
J0 = pi/4;
% g = g0*(1+0.1*randn(L,1));
% J = J0*(1+0.2*randn(L,1));
% g = g0*(2*round(rand(L,1))-1);
% J = J0*(2*round(rand(L,1))-1);
g = g0*ones(L,1);
J = J0*ones(L,1);
it = floor(L/2);

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

%% time evolution

% phi0 = rand(len,1);
% phi0 = phi0./sqrt(sum(phi0.^2));

[~,index0] = min(Hz);
phi0 = zeros(len,1);
pos = index0(1);

sz_GS0 = zeros(L,1);
temp = pos-1;
for i = 1:L
    sz_GS0(L-i+1) = mod(temp,2);
    temp = floor(temp/2);
end
sz_GS = 2*sz_GS0-1;
phi0(pos) = 1;

% flip a single point
sz_GS_new = sz_GS;
sz_GS_new(it) = - sz_GS(it);
sz_GS0_new = (sz_GS_new+1)/2;
temp = 0;
for i = 1:L
    temp = temp*2;
    temp = temp + sz_GS0_new(i);    
end
pos_new = temp+1;
phi0_e = zeros(len,1);
phi0_e(pos_new) = 1;
phi0_e = phi0_e./sqrt(sum(phi0_e.^2));

phi = phi0;
phit_store = zeros(len,nT);
phit_store(:,1) = phi0;
phi_e = phi0_e;
phit_e_store = zeros(len,nT);
phit_e_store(:,1) = phi0_e;
 
[Vx,Dx] = eig(Hx);
ex = diag(Dx);

trans_z = exp(-1i*Hz);
trans_x = exp(-1i*ex);

now = 2;
for i = 1:num_T
    % Hz
    phi = trans_z.*phi;
    phit_store(:,now) = phi;
    phi_e = trans_z.*phi_e;
    phit_e_store(:,now) = phi_e;
    now = now + 1;

    % Hx
    temp = Vx'*phi;
    trans = trans_x;
    temp = temp.*trans;
    phi = Vx*temp;
    phit_store(:,now) = phi;

    temp_e = Vx'*phi_e;
    trans = trans_x;
    temp_e = temp_e.*trans;
    phi_e = Vx*temp_e;
    phit_e_store(:,now) = phi_e;
    now = now + 1;
end

%% calculate observable

sz = zeros(L,nT);
sz_e = zeros(L,nT);
for i = 1:L
    sz(i,:) = sum(conj(phit_store).*(matrix_sz(:,i).*phit_store));
    sz_e(i,:) = sum(conj(phit_e_store).*(matrix_sz(:,i).*phit_e_store));
end
sz_diff = sz_e - sz;

sz_mean = sum(sz_GS.*sz)/L;
sz_e_mean = sum(sz_GS.*sz_e)/L;

x = (1:L)';
sz_diff_std = sqrt(sum(abs(sz_diff).*(x-it).^2));

figure;
set(gcf, 'position', [250 70 1400 900]);
subplot(1,2,1)
plot(T,sz_diff_std)

subplot(1,2,2)
% plot(log(T(1:50)),log(sz_diff_std(1:50)));
% plot(T,sz_e_mean);
imagesc(sz_diff)

toc;
