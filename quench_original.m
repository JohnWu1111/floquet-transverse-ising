clear;
clc;
format long
tic;

%% paramter
myseed = 1;
rng(myseed)

L = 4;
len = 2^L;
dt = 1e-2;
Tmax = 10;
T = 0:dt:Tmax;
nT = length(T);
g0 = 0.5;
J0 = -1;
g1 = 1;

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

[V0,D0] = eig(J0*Hxx+g0*diag(Hz)-0e-5*Hx);
phi0 = V0(:,1);
e_GS = D0(1,1);

phi = phi0;
phit_store = zeros(len,nT);
phit_store(:,1) = phi0;

[V,D] = eig(J0*Hxx+g1*diag(Hz));
e = diag(D);

trans = exp(-1i*e*dt);

for i = 2:nT
    temp = V'*phi;
    temp = temp.*trans;
    phi = V*temp;
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
plot(T,order);
% mesh(rho_xx)

toc;
