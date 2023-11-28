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
len = 2^L;
g0 = pi/2+0.0;
J0 = -1;
g1 = pi/2+0.1;
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

%% initial state
for i = 1:L/2
    Hk0 = 2*[g0+J0*cospi(k(i)) -1i*J0*sinpi(k(i));
                1i*J0*sinpi(k(i)) -g0-J0*cospi(k(i))];
    Hk0_store{i} = Hk0;
    [V,D] = eig(Hk0);
    phik0_1(i) = V(1,1);
    phik0_2(i) = V(2,1);
    e_GS(i) = D(1,1);
end
e_GS_sum = sum(e_GS);

%% time evolution
% expHk_11 = fact*(cc + 1i*ss*cospi(k));
% expHk_22 = conj(expHk_11);
% expHk_12 = -fact*ss*sinpi(k);
% expHk_21 = -conj(expHk_12);

fact = -exp(-2i*g1);
cc = cos(2*J0);
ss = sin(2*J0);

expHkz_11 = fact;
expHkz_22 = conj(expHkz_11);

expHkxx_11 = cc - 1i*ss*cospi(k);
expHkxx_22 = conj(expHkxx_11);
expHkxx_12 = -ss*sinpi(k);
expHkxx_21 = -conj(expHkxx_12);
norm_expHkxx = expHkxx_11.*expHkxx_22 - expHkxx_12.*expHkxx_21;

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
    phik_1t = expHkxx_11.*phik_1 + expHkxx_12.*phik_2;
    phik_2 = expHkxx_21.*phik_1 + expHkxx_22.*phik_2;
    phik_1 = phik_1t;
    phit_1_store(:,2*i) = phik_1;
    phit_2_store(:,2*i) = phik_2;

    n(2*i) = 2*(phik_1'*phik_1)/L;
    Cr(:,2*i) = 2*(abs(phik_1).^2)'*cospi(k*x')/L;
    Fr(:,2*i) = -2i*(phik_1.*conj(phik_2))'*sinpi(k*x')/L;
    n1n2(2*i) = n(2*i)^2 + Fr(2)*conj(Fr(2)) - Cr(2)^2;
    sz1sz2(2*i) = 1-4*n(2*i)+4*n1n2(2*i);

    phik_1 = expHkz_11.*phik_1;
    phik_2 = expHkz_22.*phik_2;
    phit_1_store(:,2*i+1) = phik_1;
    phit_2_store(:,2*i+1) = phik_2;

    n(2*i+1) = 2*(phik_1'*phik_1)/L;
    Cr(:,2*i+1) = 2*(abs(phik_1).^2)'*cospi(k*x')/L;
    Fr(:,2*i+1) = -2i*(phik_1.*conj(phik_2))'*sinpi(k*x')/L;
    n1n2(2*i+1) = n(2*i+1)^2 + Fr(2)*conj(Fr(2)) - Cr(2)^2;
    sz1sz2(2*i+1) = 1-4*n(2*i+1)+4*n1n2(2*i+1);
end

%% calculate observable

figure;
set(gcf, 'position', [250 70 1400 900]);
% subplot(2,1,1)
% plot(T,sz1sz2);
plot(T,n);

toc;