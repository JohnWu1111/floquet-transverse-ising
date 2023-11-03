clear;
clc;
format long
tic;

%% paramter
myseed = 1;
rng(myseed)

L = 10;
num_T = 200;
T = 0:2*num_T;
nT = 2*num_T+1;
dt = 1;
len = 2^L;
g = 1;
J = -1;
g0 = 0.01;

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

[V0,D0] = eig(J*Hxx+g0*diag(Hz));
phi0 = V0(:,1);
e_GS = D0(1,1);

phi = phi0;
phit_store = zeros(len,nT);
phit_store(:,1) = phi0;

[V1,D1] = eig(J*Hxx+g*diag(Hz));
e1 = diag(D1);
[V2,D2] = eig(J*Hxx-g*diag(Hz));
e2 = diag(D2);

trans1 = exp(-1i*e1);
trans2 = exp(-1i*e2);

for i = 1:num_T
    % Hxx+Hz
    temp = V1'*phi;
    trans = trans1;
    temp = temp.*trans;
    phi = V1*temp;
    phit_store(:,2*i) = phi;

    % Hxx-Hz
    temp = V2'*phi;
    trans = trans2;
    temp = temp.*trans;
    phi = V2*temp;
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
