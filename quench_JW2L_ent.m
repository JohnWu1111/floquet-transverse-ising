% clear;
clc;
format long
tic;

%% paramter
myseed = 1;
rng(myseed)

L = 6;
dt = 1e-2;
Tmax = 10;
T = 0:dt:Tmax;
nT = length(T);
g0 = 0;
J0 = -1;
g1 = 1;
x = (0:L)';
xx = (1:L)';

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
Green0 = V_GS*V_GS';
% V_edge = [V(:,1:L-1),V(:,L+1)];
% G0 = V_edge*V_edge';
Green = Green0;
Greent_store = cell(nT,1);
Greent_store{1} = Green;
% DagDag = V_GS*circshift(V_GS,L,1)';

M_fact = cell(L-1,1);
sigma_y = sparse([0 -1i;1i 0]);
M_fact{1} = sigma_y.';
for i = 1:L-1
    In = sparse(eye(i));
    M_fact{i} = kron(sigma_y,In).';
end

%% time evolution
% Hxx = [J0*A,J0*B;
%       -J0*B,-J0*A];
[V,D] = eig(g1*diag(Hz) + J0*Hxx);
e = diag(D);

trans = exp(2i*e*dt);
trans0 = V.*trans'*V';

sz1sz2 = zeros(1,nT);
sz1sz2(1) = 1-4*Green(1,1)+4*(Green(1,1)^2+Green(1+L,2)^2-Green(1,2)^2);
n = zeros(L,nT);
n0 = diag(Green);
n(:,1) = n0(1:L);
ent = zeros(1,nT);
ent(1) = cal_ent2(Green,L);
% ent(1) = 0;


for i = 2:nT
    Green = trans0'*Green*trans0;
    Greent_store{i} = Green;
    sz1sz2(i) = real(1-4*Green(1,1)+4*(Green(1,1)^2+abs(Green(1+L,2))^2-Green(1,2)^2));
%     n0 = real(diag(Green));
%     n(:,i) = n0(1:L);
    ent(i) = cal_ent2(Green,L);
end


%% calculate observable
figure;
set(gcf, 'position', [250 70 1400 900]);
% plot(xx,n(:,1));
% plot(T,n(100,:))
plot(T,real(ent))
% mesh(rho_xx)

toc;


% function ent = cal_ent(Green,L)
% Green_r = Green(1:L,1:L);
% DagDag_r = Green(1:L,L+1:2*L);
% 
% x = 1:L/2;
% x = 2*x-1;
% C_odd = Green_r(x,x);
% F_odd = DagDag_r(x,x);
% 
% ent_M = [eye(L/2)-C_odd,F_odd;
%          F_odd',C_odd];
% e = 2*eig(ent_M);
% ent = -sum(e.*log(e));
% end

function ent = cal_ent2(Green,L)
Green_r = Green(1:L,1:L);
DagDag_r = Green(1:L,L+1:2*L);

x = 1:L/2;
% x = 2*x;
C_odd = Green_r(x,x);
F_odd = DagDag_r(x,x);

ent_M = (2*C_odd-eye(L/2)-2*F_odd)*(2*C_odd-eye(L/2)+2*conj(F_odd));
e = eig(ent_M);
ep = log(2*sqrt(1./(1-e))+2*sqrt(1./(1-e)-4));
rho = -ep./sum(ep);
ent = -sum(rho.*exp(rho));
end
