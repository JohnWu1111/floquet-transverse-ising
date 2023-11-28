clear;
clc;
format long
tic;

%% paramter

L = 1000;
num_T = 100;
T = 0:num_T;
nT = num_T+1;
len = 2^L;
g0 = 0;
J = -1;
g = 1;
% dt = 1*pi/(2*sqrt(J^2+g^2));
% dt = 2*pi/(2*J);
dt = 3.1;
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
n = zeros(L,nT);
n0 = diag(G);
n(:,1) = n0(1:L);
[Vp,Dp] = eig(J*Hxx+g*diag(Hz));
ep = diag(Dp);
[Vm,Dm] = eig(J*Hxx-g*diag(Hz));
em = diag(Dm);

trans_p = exp(2i*ep*dt);
trans0_p = Vp.*trans_p'*Vp';
trans_m = exp(2i*em*dt);
trans0_m = Vm.*trans_m'*Vm';
trans = trans0_p*trans0_m;
H_eff = -1i*logm(trans);

[V_eff,D_eff] = eig(H_eff);
e_eff = real(diag(D_eff));
[e_eff1,I] = sort(e_eff);
% e_eff = e_eff(I);
V_eff = V_eff(:,I);
e_eff = diag(D_eff);
VeffGS = V_eff(:,1:L);
% VeffGS = V_eff(:,[1:L-1,L+1]);
Geff0 = VeffGS*VeffGS';
neff0 = diag(Geff0);
neff = neff0(1:L);

figure;
set(gcf, 'position', [250 70 1400 900]);
% plot(T,sz1sz2);
% mesh(n)
plot(1:L,neff)

for i = 2:nT
    G = trans'*G*trans;
    Gt_store{2*i} = G;
    sz1sz2(2*i) = 1-4*G(1,1)+4*(G(1,1)^2+abs(G(1+L,2))^2-G(1,2)^2);
    n0 = real(diag(G));
    n(:,i) = n0(1:L);
end

%% calculate observable

figure;
set(gcf, 'position', [250 70 1400 900]);
% plot(T,sz1sz2);
% mesh(n)
plot(1:L,neff)

toc;