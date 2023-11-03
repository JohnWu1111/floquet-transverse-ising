clear;
clc;
format long
tic;

%% paramter
myseed = 1;
rng(myseed)

L = 8;
num_T = 200;
T = 0:2*num_T;
nT = 2*num_T+1;
len = 2^L;
g0 = pi/2+0.0;
J0 = -1;
g1 = pi/2+0.1;

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
matrix_sx = cell(L,1);
for i = 1:L
    Hz_temp = ones(2^(i-1),1);
    Hz_temp = kron(Hz_temp,sigmaz);
    Hz_temp = kron(Hz_temp,ones(2^(L-i),1));
    Hz = Hz + Hz_temp;

    sx_temp = eye(2^(i-1));
    sx_temp = kron(sx_temp,sigmax);
    sx_temp = kron(sx_temp,eye(2^(L-i)));
    matrix_sx{i} = sx_temp;
end

%% time evolution

[V0,D0] = eig(J0*Hxx+g0*diag(Hz));
phi0 = V0(:,1);
e_GS = D0(1,1);

phi = phi0;
phit_store = zeros(len,nT);
phit_store(:,1) = phi0;

[Vx,Dx] = eig(J0*Hxx);
ex = diag(Dx);

trans_z = exp(-1i*g1*Hz);
trans_x = exp(-1i*ex);

for i = 1:num_T
    % Hxx
    temp = Vx'*phi;
    trans = trans_x;
    temp = temp.*trans;
    phi = Vx*temp;
    phit_store(:,2*i) = phi;

    % Hz
    phi = trans_z.*phi;
    phit_store(:,2*i+1) = phi;
end

%% calculate observable

Msz1sz2 = kron(sigmaz,sigmaz);
Msz1sz2 = kron(Msz1sz2,ones(2^(L-2),1));
sz1sz2 = sum(conj(phit_store).*(Msz1sz2.*phit_store));

% sx = zeros(L,nT);
% for i = 1:L
%     sx(i,:) = sum(conj(phit_store).*(matrix_sx{i}*phit_store));
% end
% sx_mean = sum(sx)/L;

figure;
set(gcf, 'position', [250 70 1400 900]);
% subplot(2,1,1)
plot(T,sz1sz2);
% subplot(2,1,2)
% plot(T,sx_mean);

toc;
