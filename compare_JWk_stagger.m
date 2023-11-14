clear;
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
k = (1/L:2/L:(L-1)/L)';
k_full = (-(L-1)/L:2/L:(L-1)/L)';
k0 = (0:2/L:(L-2)/L)';
x = (0:L-1)';
xx = (1:L)';

% expHk_store = cell(L/2,1);
Hk0_store = cell(L/2,1);
phik0_1 = zeros(L/2,1);
phik0_2 = zeros(L/2,1);
e_GS = zeros(L/2,1);

ck = cospi(k);
sk = sinpi(k);

%% initial state
for i = 1:L/2
    Hk0 = 2*[g0+J*ck(i) -1i*J*sk(i);
                1i*J*sk(i) -g0-J*ck(i)];
    Hk0_store{i} = Hk0;
    [V,D] = eig(Hk0);
    phik0_1(i) = V(1,1);
    phik0_2(i) = V(2,1);
    e_GS(i) = D(1,1);
end
e_GS_sum = sum(e_GS);

%% time evolution

fact = sqrt(g^2+J^2+2*g*J*ck);
sf = sin(2*fact*dt)./fact;

expHkp_11 = cos(2*dt*fact) - 1i*(g+J*ck).*sf;
expHkp_22 = conj(expHkp_11);
expHkp_12 = -J*sk.*sf;
expHkp_21 = -expHkp_12;
norm_expHkp = abs(expHkp_11.*expHkp_22 - expHkp_12.*expHkp_21);

fact = sqrt(g^2+J^2-2*g*J*ck);
sf = sin(2*fact*dt)./fact;

expHkm_11 = cos(2*dt*fact) - 1i*(-g+J*ck).*sf;
expHkm_22 = conj(expHkm_11);
expHkm_12 = -J*sk.*sf;
expHkm_21 = -expHkm_12;
norm_expHkm = abs(expHkm_11.*expHkm_22 - expHkm_12.*expHkm_21);

phik_1 = phik0_1;
phik_2 = phik0_2;
phit_1_store = zeros(L/2,nT);
phit_2_store = zeros(L/2,nT);
phit_1_store(:,1) = phik0_1;
phit_2_store(:,1) = phik0_2;

sz1sz2 = zeros(1,nT);
Cr = zeros(L,nT);
Fr = zeros(L,nT);
n = zeros(1,nT);
n1n2 = zeros(1,nT);
n(1) = 2*(phik_1'*phik_1)/L;
Cr(:,1) = 2*(abs(phik_1).^2)'*cospi(k*x')/L;
Fr(:,1) = -2i*(phik_1.*conj(phik_2))'*sinpi(k*x')/L;
n1n2(1) = n(1)^2 + Fr(2)^2 - Cr(2)^2;
sz1sz2(1) = 1-4*n(1)+4*n1n2(1);

for i = 1:num_T
    phik_1t = expHkp_11.*phik_1 + expHkp_12.*phik_2;
    phik_2 = expHkp_21.*phik_1 + expHkp_22.*phik_2;
    phik_1 = phik_1t;
    phit_1_store(:,2*i) = phik_1;
    phit_2_store(:,2*i) = phik_2;

    n(2*i) = 2*(phik_1'*phik_1)/L;
    Cr(:,2*i) = 2*(abs(phik_1).^2)'*cospi(k*x')/L;
    Fr(:,2*i) = -2i*(phik_1.*conj(phik_2))'*sinpi(k*x')/L;
    n1n2(2*i) = n(2*i)^2 + Fr(2,2*i)*conj(Fr(2,2*i)) - Cr(2,2*i)^2;
    sz1sz2(2*i) = 1-4*n(2*i)+4*n1n2(2*i);

    phik_1t = expHkm_11.*phik_1 + expHkm_12.*phik_2;
    phik_2 = expHkm_21.*phik_1 + expHkm_22.*phik_2;
    phik_1 = phik_1t;
    phit_1_store(:,2*i) = phik_1;
    phit_2_store(:,2*i) = phik_2;

    n(2*i+1) = 2*(phik_1'*phik_1)/L;
    Cr(:,2*i+1) = 2*(abs(phik_1).^2)'*cospi(k*x')/L;
    Fr(:,2*i+1) = -2i*(phik_1.*conj(phik_2))'*sinpi(k*x')/L;
    n1n2(2*i+1) = n(2*i+1)^2 + Fr(2,2*i+1)*conj(Fr(2,2*i+1)) - Cr(2,2*i+1)^2;
    sz1sz2(2*i+1) = 1-4*n(2*i+1)+4*n1n2(2*i+1);
end

%% calculate observable

figure;
set(gcf, 'position', [250 70 1400 900]);
% subplot(2,1,1)
plot(T,sz1sz2);

toc;

% expHk_11 = fact*(cc + 1i*ss*cospi(k));
% expHk_22 = conj(expHk_11);
% expHk_12 = -fact*ss*sinpi(k);
% expHk_21 = -conj(expHk_12);

% expHkp_11 = zeros(L/2,1);
% expHkp_22 = zeros(L/2,1);
% expHkp_12 = zeros(L/2,1);
% expHkp_21 = zeros(L/2,1);
% expHkm_11 = zeros(L/2,1);
% expHkm_22 = zeros(L/2,1);
% expHkm_12 = zeros(L/2,1);
% expHkm_21 = zeros(L/2,1);
% expHkp_store = cell(L/2,1);
% expHkm_store = cell(L/2,1);
% for i = 1:L/2
%     Hkp = 2*[g+J*ck(i) -1i*J*sk(i);
%                 1i*J*sk(i) -g-J*ck(i)];
%     Hkm = 2*[-g+J*ck(i) -1i*J*sk(i);
%         1i*J*sk(i) g-J*ck(i)];
%     expHkp = expm(-1i*Hkp*dt);
%     expHkm = expm(-1i*Hkm*dt);
%     expHkp_store{i} = expHkp;
%     expHkm_store{i} = expHkm;
%     expHkp_11(i) = expHkp(1,1);
%     expHkp_22(i) = expHkp(2,2);
%     expHkp_12(i) = expHkp(1,2);
%     expHkp_21(i) = expHkp(2,1);
%     expHkm_11(i) = expHkm(1,1);
%     expHkm_22(i) = expHkm(2,2);
%     expHkm_12(i) = expHkm(1,2);
%     expHkm_21(i) = expHkm(2,1);
% end