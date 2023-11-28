% clear;
clc;
format long
tic;

%% paramter
myseed = 1;
rng(myseed)

L = 50;
dt = 20;
Tmax = 20000;
T = 0:dt:Tmax;
nT = length(T);
g0 = 0.0;
J0 = -1;
g1 = 1;
g2 = -1;
dt0 = 0.5*dt;
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

rho_xx = zeros(L,nT);
rho_xx(:,1) = cal_rho_xx(Green,M_fact,L);

%% time evolution
% Hxx = [J0*A,J0*B;
%       -J0*B,-J0*A];
[V1,D1] = eig(g1*diag(Hz) + J0*Hxx);
e1 = diag(D1);
[V2,D2] = eig(g2*diag(Hz) + J0*Hxx);
e2 = diag(D2);

trans1 = exp(2i*e1*(dt-dt0));
trans0_1 = V1.*trans1'*V1';
trans2 = exp(2i*e2*dt);
trans0_2 = V2.*trans2'*V2';
trans3 = exp(2i*e1*dt0);
trans0_3 = V1.*trans2'*V1';

% sz1sz2 = zeros(1,nT);
% sz1sz2(1) = 1-4*Green(1,1)+4*(Green(1,1)^2+Green(1+L,2)^2-Green(1,2)^2);
% n = zeros(L,nT);
% n0 = diag(Green);
% n(:,1) = n0(1:L);

for i = 2:nT
    Green = trans0_1'*Green*trans0_1;
    Green = trans0_2'*Green*trans0_2;
    Green = trans0_3'*Green*trans0_3;
    Greent_store{i} = Green;
%     sz1sz2(i) = real(1-4*Green(1,1)+4*(Green(1,1)^2+abs(Green(1+L,2))^2-Green(1,2)^2));
%     n0 = real(diag(Green));
%     n(:,i) = n0(1:L);
    rho_xx(:,i) = cal_rho_xx(Green,M_fact,L);
end

order = sqrt(sum(rho_xx)/L);

%% calculate observable
figure;
set(gcf, 'position', [250 70 1400 900]);
% plot(xx,n(:,1));
% plot(T,n(100,:))
% plot(T,order)
mesh(rho_xx)

toc;


function rho_xx = cal_rho_xx(Green,M_fact,L)
Green_r = Green(1:L,1:L);
DagDag_r = Green(1:L,L+1:2*L);

Green_r = circshift(Green_r,1,1);
DagDag_r = circshift(DagDag_r,1,1);
Green_r(1,:) = -Green_r(1,:);
DagDag_r(1,:) = -DagDag_r(1,:);

S = 2i*imag(DagDag_r);
G = 2*real(Green_r + DagDag_r);
G = G -circshift(eye(L),1);

S = triu(circshift(S,-1),1);
% S = S - triu(S,1)';

rho_xx = zeros(L,1);
rho_xx(1) = 1;

for i = 2:L
    n = i-1;
    range = 1:n;
    PF = [S(range,range), G(range,range);
          zeros(n),S(range,range)];
    PF = PF - triu(PF,1).';
    e = eig(M_fact{n}*PF);
    rho_xx(i) = real(1i^(n)*exp(sum(log(e))/2));
%     rho_xx(i) = pfaffian_LTL(PF);
%     rho_xx(i) = 1i^(n^2)*exp(trace(logm(M_fact{n}*PF))/2);
%     rho_xx(i) = det(G(range,range));
end
end
