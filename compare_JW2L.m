clear;
clc;
format long
tic;

%% paramter
myseed = 1;
rng(myseed)

L = 200;
num_T = 200;
T = 0:2*num_T;
nT = 2*num_T+1;
g0 = pi/2+0.0;
J0 = -1;
g1 = pi/2+0.1;
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
[V, D] = eig(H0);
e = diag(D);
e_GS = sum(e(1:L));

V_GS = V(:,1:L);
G0 = V_GS*V_GS';
G = G0;
Gt_store = cell(nT,1);
Gt_store{1} = G;

%% time evolution
% Hxx = [J0*A,J0*B;
%       -J0*B,-J0*A];
[Vx,Dx] = eig(J0*Hxx);
ex = diag(Dx);

% Hz = [g1*ones(L,1);-g1*ones(L,1)];
trans_z = exp(2i*g1*Hz);
trans_x = exp(2i*ex);
% trans0_x = (V*diag(trans_x))/V;
% trans0_x = expm(-1i*J0*Hxx*2);
trans0_x = Vx.*trans_x'*Vx';
% trans0_x = Vx*diag(trans_x)'*Vx';

sz1sz2 = zeros(1,nT);
sz1sz2(1) = 1-4*G(1,1)+4*(G(1,1)^2+G(1+L,2)^2-G(1,2)^2);
n1 = zeros(1,nT);
n1(1) = G(1,1);

for i = 1:num_T
    % Hxx
%     GG = Vx'*G*Vx;
%     G_temp = trans_x'.*GG.*trans_x;
%     G = Vx*G_temp*Vx';

    G = trans0_x'*G*trans0_x;

    Gt_store{2*i} = G;
    sz1sz2(2*i) = 1-4*G(1,1)+4*(G(1,1)^2+abs(G(1+L,2))^2-G(1,2)^2);
    n1(2*i) = G(1,1);

    % Hz
    G = trans_z'.*G.*trans_z;
    Gt_store{2*i+1} = G;
    sz1sz2(2*i+1) = 1-4*G(1,1)+4*(G(1,1)^2+abs(G(1+L,2))^2-G(1,2)^2);
    n1(2*i+1) = G(1,1);
end

%% calculate observable
figure;
set(gcf, 'position', [250 70 1400 900]);
% subplot(2,1,1)
plot(T,n1);

toc;