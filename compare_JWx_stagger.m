% clear;
clc;
format long
tic;

%% paramter
myseed = 1;
rng(myseed)

L = 4;
num_T = 200;
T = 0:2*num_T;
nT = 2*num_T+1;
dt = 1;
len = 2^L;
g0 = 0.1;
J = -1;
g = 1;
x = (0:L)';
xx = (1:L)';

C_dag = [0 1;0 0];
C = [0 0;1 0];
ni = [1;0];
K = [1 0;0 -1];
I2 = eye(2);
sigmaz = [1;-1];

%% construction of Hamiltonian
% Hxx
Hxx = zeros(len);
for i = 1:L-1
    temp = 1;
    temp = kron(temp,eye(2^(i-1)));
    temp = kron(temp,(C_dag+C));
    temp = kron(temp,(C_dag+C));
    temp = kron(temp,eye(2^(L-i-1)));
    Hxx = Hxx + temp;
end
% i = L
temp = C_dag-C;
for i = 2:L-1
    temp = kron(temp,K);
end
temp = kron(temp,(C_dag-C));
Hxx = Hxx + temp;

% Hz
Hz = zeros(len,1);
for i = 1:L
    temp = 1;
    temp = kron(temp, ones(2^(i-1),1));
    temp = kron(temp,ni);
    temp = kron(temp, ones(2^(L-i),1));
    Hz = Hz + temp;
end
Hz = 2*Hz - L;

%% initial state

H0 = g0*diag(Hz) + J*Hxx;
[V, D] = eig(H0);
phi0 = V(:,1);
e_GS = D(1,1);

phi = phi0;
phit_store = zeros(len,nT);
phit_store(:,1) = phi0;

H1 = J*Hxx+g*diag(Hz);
[V1,D1] = eig(H1);
e1 = diag(D1);
H2 = J*Hxx-g*diag(Hz);
[V2,D2] = eig(H2);
e2 = diag(D2);

trans1 = exp(-1i*e1*dt);
trans2 = exp(-1i*e2*dt);

for i = 1:num_T
    % Hxx+Hz
    temp = V1'*phi;
    temp = temp.*trans1;
    phi = V1*temp;
    phit_store(:,2*i) = phi;

    % Hxx-Hz
    temp = V2'*phi;
    temp = temp.*trans2;
    phi = V2*temp;
    phit_store(:,2*i+1) = phi;
end

%% calculate observable

Msz1sz2 = kron(sigmaz,sigmaz);
Msz1sz2 = kron(Msz1sz2,ones(2^(L-2),1));
Mn1 = kron(ni,ones(2^(L-1),1));
Mn2 = kron([1;1],ni);
Mn2 = kron(Mn2,ones(2^(L-2),1));
Mn3 = kron(ones(4,1),ni);
Mn3 = kron(Mn3,ones(2^(L-3),1));
Mn1n2 = kron(ni,ni);
Mn1n2 = kron(Mn1n2,ones(2^(L-2),1));
sz1sz2 = sum(conj(phit_store).*(Msz1sz2.*phit_store));
n1 = sum(conj(phit_store).*(Mn1.*phit_store));
n2 = sum(conj(phit_store).*(Mn2.*phit_store));
n3 = sum(conj(phit_store).*(Mn3.*phit_store));
n1n2 = sum(conj(phit_store).*(Mn1n2.*phit_store));

MC1j = cell(L,1);
MC1j{1} = diag(Mn1);
C1j = zeros(L,nT);
C1j(1,:) = n1;
MF1j = cell(L,1);
F1j = zeros(L,nT);
for j = 2:L
    temp = C_dag;
    for i = 1:j-2
        temp = kron(temp,K);
    end
%     temp = kron(temp,eye(2^(j-2)));
    temp = kron(temp,C);
    temp = kron(temp,eye(2^(L-j)));
    MC1j{j} = temp;
    C1j(j,:) = sum(conj(phit_store).*(temp*phit_store));

    temp = C_dag;
    for i = 1:j-2
        temp = kron(temp,K);
    end
%     temp = kron(temp,eye(2^(j-2)));
    temp = kron(temp,C_dag);
%     temp = kron(temp,eye(2^(L-j)));
    for i = j+1:L
        temp = kron(temp,K);
    end
    MF1j{j} = temp;
    F1j(j,:) = sum(conj(phit_store).*(temp*phit_store));
end

n1n2_check = n1.*n2+abs(F1j(2,:)).^2-abs(C1j(2,:)).^2;

% sx = zeros(L,nT);
% for i = 1:L
%     sx(i,:) = sum(conj(phit_store).*(matrix_sx{i}*phit_store));
% end
% sx_mean = sum(sx)/L;

figure;
set(gcf, 'position', [250 70 1400 900]);
% subplot(2,1,1)
plot(T,sz1sz2);

toc;