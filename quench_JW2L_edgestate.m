clear;
clc;
format long
tic;

%% paramter
myseed = 1;
rng(myseed)

L = 100;
dt = 1e-2;
Tmax = 100;
T = 0:dt:Tmax;
nT = length(T);
g0 = 0;
J0 = -1;
g1 = 0.5;
x = (0:L)';
xx = (1:L)';

%% construction of Hamiltonian
A = zeros(L);
for i = 1:L-1
    A(i,i+1) = 1/2;
    A(i+1,i) = 1/2;
end
% A(L,1) = -1/2;
% A(1,L) = -1/2;

B = zeros(L);
for i = 1:L-1
    B(i,i+1) = 1/2;
    B(i+1,i) = -1/2;
end
% B(L,1) = -1/2;
% B(1,L) = 1/2;

%% initial state
% H0 = [g0*eye(L)+J0*A,J0*B;
%       -J0*B,-g0*eye(L)-J0*A];
Hxx = [A,B;
      -B,-A];
Hz = [ones(L,1);-ones(L,1)];
H0 = g0*diag(Hz) + J0*Hxx;
[V0, D0] = eig(H0);
e0 = diag(D0);
e_GS = sum(e0(1:L));

V_GS = V0(:,1:L);
G0 = V_GS*V_GS';
% V_edge = [V(:,1:L-1),V(:,L+1)];
% G0 = V_edge*V_edge';
G = G0;
Gt_store = cell(nT,1);
Gt_store{1} = G;

%% time evolution
% Hxx = [J0*A,J0*B;
%       -J0*B,-J0*A];
[V,D] = eig(g1*diag(Hz) + J0*Hxx);
e = diag(D);

trans = exp(2i*e*dt);
trans0 = V.*trans'*V';

sz1sz2 = zeros(1,nT);
sz1sz2(1) = 1-4*G(1,1)+4*(G(1,1)^2+G(1+L,2)^2-G(1,2)^2);
n = zeros(L,nT);
n0 = diag(G);
n(:,1) = n0(1:L);

for i = 2:nT
    G = trans0'*G*trans0;
    Gt_store{2*i} = G;
    sz1sz2(i) = 1-4*G(1,1)+4*(G(1,1)^2+abs(G(1+L,2))^2-G(1,2)^2);
    n0 = real(diag(G));
    n(:,i) = n0(1:L);
end

%% calculate observable
figure;
set(gcf, 'position', [250 70 1400 900]);
% plot(xx,n(:,1));
% plot(T,n(100,:))
mesh(n)

toc;