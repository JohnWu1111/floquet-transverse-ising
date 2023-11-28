clear;
clc;
format long
tic;

%% paramter
myseed = 1;
rng(myseed)

L = 6;
Tmax = 100;
dt = 1;
T = 0:dt:Tmax;
nT = length(T);

len = 2^L;
g = 1;
J = -1;
g0 = 0.5;

sigmaz = [1;-1];
sigmax = [0 1;1 0];
I2 = eye(2);
II2 = ones(2);

%% construction of Hamiltonian and observable
Hxx = zeros(len);
for i = 1:L-1
    Hxx_temp = eye(2^(i-1));
    Hxx_temp = kron(Hxx_temp, sigmax);
    Hxx_temp = kron(Hxx_temp, sigmax);
    Hxx_temp = kron(Hxx_temp, eye(2^(L-i-1)));
    Hxx = Hxx + Hxx_temp;
end
% PBC
Hxx_temp = sigmax;
Hxx_temp = kron(Hxx_temp, eye(2^(L-2)));
Hxx_temp = kron(Hxx_temp, sigmax);
Hxx = Hxx + Hxx_temp;

Hz = zeros(len,1);
Hx = zeros(len);
matrix_sx = cell(L,1);
matrix_sz = zeros(len,L);
for i = 1:L
    Hz_temp = ones(2^(i-1),1);
    Hz_temp = kron(Hz_temp,sigmaz);
    Hz_temp = kron(Hz_temp,ones(2^(L-i),1));
    Hz = Hz + Hz_temp;
    Hx_temp = eye(2^(i-1));
    Hx_temp = kron(Hx_temp,sigmax);
    Hx_temp = kron(Hx_temp,eye(2^(L-i)));
    Hx = Hx + Hx_temp;

    sx_temp = eye(2^(i-1));
    sx_temp = kron(sx_temp,sigmax);
    sx_temp = kron(sx_temp,eye(2^(L-i)));
    matrix_sx{i} = sx_temp;

    sz_temp = ones(2^(i-1),1);
    sz_temp = kron(sz_temp,sigmaz);
    sz_temp = kron(sz_temp,ones(2^(L-i),1));
    matrix_sz(:,i) = sz_temp;
end

%% time evolution

[V0,D0] = eig(J*Hxx+g0*diag(Hz));
phi0 = V0(:,1);
e_GS = D0(1,1);

phi = phi0;
phit_store = zeros(len,nT);
phit_store(:,1) = phi0;

H1 = J*Hxx+g*diag(Hz);
[V1,D1] = eig(H1);
e1 = diag(D1);
[V2,D2] = eig(J*Hxx-g*diag(Hz));
e2 = diag(D2);

trans1 = exp(-1i*e1*dt);
trans2 = exp(-1i*e2*dt);

for i = 2:nT
    % Hxx+Hz
    temp = V1'*phi;
    temp = temp.*trans1;
    phi = V1*temp;

    % Hxx-Hz
    temp = V2'*phi;
    temp = temp.*trans2;
    phi = V2*temp;
    phit_store(:,i) = phi;
end

%% calculate observable

Msz1sz2 = kron(sigmaz,sigmaz);
Msz1sz2 = kron(Msz1sz2,ones(2^(L-2),1));
sz1sz2 = sum(conj(phit_store).*(Msz1sz2.*phit_store));

rho_xx = zeros(L,nT);
rho_xx(1,:) = 1;
for i = 2:L
    Msxsx = sigmax;
    Msxsx = kron(Msxsx,eye(2^(i-2)));
    Msxsx = kron(Msxsx,sigmax);
    Msxsx = kron(Msxsx,eye(2^(L-i)));
    rho_xx(i,:) = real(sum(conj(phit_store).*(Msxsx*phit_store)));
end

sx = zeros(L,nT);
sz = zeros(L,nT);
for i = 1:L
    sx(i,:) = real(sum(conj(phit_store).*(matrix_sx{i}*phit_store)));
    sz(i,:) = sum(conj(phit_store).*matrix_sz(:,i).*phit_store);
%     for j = 1:nT
%         sx(i,j) = phit_store(:,j)'*matrix_sx{i}*phit_store(:,j);
%     end
end
order = sqrt(sum(rho_xx)/L);

figure;
set(gcf, 'position', [250 70 1400 900]);
% subplot(2,1,1)
% plot(T,order);
mesh(rho_xx)
% subplot(2,1,2)
% plot(T,sx_mean);

toc;
