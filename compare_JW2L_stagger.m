clear;
clc;
format long
tic;

%% paramter

L = 4;
num_T = 200;
T = 0:2*num_T;
nT = 2*num_T+1;
dt = 1;
T_eff = 0:num_T;
nT_eff = num_T+1;
len = 2^L;
g0 = 0.1;
J = -1;
g = 1;
% k = (2/L:2/L:1)';
k = (1/L:2/L:(L-1)/L)';
x = (0:L)';
xx = (0:L-1)';

%% construction of Hamiltonian
A = zeros(L);
for i = 1:L-1
    A(i,i+1) = 1/2;
    A(i+1,i) = 1/2;
end
A(L,1) = -1/2;
A(1,L) = -1/2;

B = zeros(L);
for i = 1:L-1
    B(i,i+1) = 1/2;
    B(i+1,i) = -1/2;
end
B(L,1) = -1/2;
B(1,L) = 1/2;

%% initial state
Hxx = [A,B;
      -B,-A];
Hz = [ones(L,1);-ones(L,1)];
H0 = g0*diag(Hz) + J*Hxx;
[V, D] = eig(H0);
e = diag(D);
e_GS = sum(e(1:L));

V_GS = V(:,1:L);
G0 = V_GS*V_GS';
G = G0;
Gt_store = cell(nT,1);
Gt_store{1} = G;

%% time evolution
sz1sz2 = zeros(1,nT);
sz1sz2(1) = 1-4*G(1,1)+4*(G(1,1)^2+G(1+L,2)^2-G(1,2)^2);
n1 = zeros(1,nT);
n1(1) = G(1,1);
[Vp,Dp] = eig(J*Hxx+g*diag(Hz));
ep = diag(Dp);
[Vm,Dm] = eig(J*Hxx-g*diag(Hz));
em = diag(Dm);

trans_p = exp(2i*ep);
trans0_p = Vp.*trans_p'*Vp';
trans_m = exp(2i*em);
trans0_m = Vm.*trans_m'*Vm';

for i = 1:num_T
    % Hp
    G = trans0_p'*G*trans0_p;
    Gt_store{2*i} = G;
    sz1sz2(2*i) = 1-4*G(1,1)+4*(G(1,1)^2+abs(G(1+L,2))^2-G(1,2)^2);
    n1(2*i) = G(1,1);

    % Hm
    G = trans0_m'*G*trans0_m;
    Gt_store{2*i+1} = G;
    sz1sz2(2*i+1) = 1-4*G(1,1)+4*(G(1,1)^2+abs(G(1+L,2))^2-G(1,2)^2);
    n1(2*i+1) = G(1,1);
end

%% calculate observable

figure;
set(gcf, 'position', [250 70 1400 900]);
% subplot(2,1,1)
plot(T,sz1sz2);

toc;