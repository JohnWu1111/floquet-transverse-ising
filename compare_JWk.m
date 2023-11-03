clear;
clc;
format long
tic;

%% paramter
myseed = 1;
rng(myseed)

L = 8;
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

fact = exp(-2i*(g1-pi/2));
cc = cos(2*J0);
ss = sin(2*J0);
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

expHkz_11 = fact;
expHkz_22 = conj(expHkz_11);

expHkxx_11 = cc + 1i*ss*cospi(k);
expHkxx_22 = conj(expHkxx_11);
expHkxx_12 = -ss*sinpi(k);
expHkxx_21 = -conj(expHkxx_12);

phik_1 = phik0_1;
phik_2 = phik0_2;
phit_1_store = zeros(L/2,nT);
phit_2_store = zeros(L/2,nT);
phit_1_store(:,1) = phik0_1;
phit_2_store(:,1) = phik0_2;

sz1sz2_t = zeros(1,nT);
n = 2*(phik_1'*phik_1)/L;
Cr = 2*(abs(phik_1).^2)'*cospi(k*x')/L;
Fr = -2i*(phik_1.*conj(phik_2))'*sinpi(k*x')/L;
% % n2 = 2*(phik_1'*(cospi(k-k')+cospi(k+k'))*phik_1)/(L^2);
% temp1 = sum((phik_1.*cospi(k))'*phik_2);
% temp2 = sum(phik_2'*(cospi(k).*phik_1));
% temp3 = 2*sum(cospi(k).*abs(phik_1).^2)/L;
% temp4 = 2*sum(cospi(k).*abs(phik_2).^2)/L;
% n2 = n^2 - temp1*temp2 + temp3*temp4;


sz1sz2_t(1) = 1-4*n+4*n2;

for i = 1:num_T
    phik_1t = expHkxx_11.*phik_1 + expHkxx_12.*phik_2;
    phik_2 = expHkxx_21.*phik_1 + expHkxx_22.*phik_2;
    phik_1 = phik_1t;
    phit_1_store(:,2*i) = phik_1;
    phit_2_store(:,2*i) = phik_2;

    phik_full = [phik_1;flip(conj(phik_2))];
    temp1 = exp(-1i*k_full*pi).*phik_full;
    temp2 = exp(-2i*k_full*pi).*phik_full;
    n1 = sum(temp1*temp1',"all")/L;
    n2 = sum(temp2*temp2',"all")/L;
    sz1sz2_t(2*i+1) = (1-2*n1)*(1-2*n2);

    phik_1 = expHkz_11.*phik_1;
    phik_2 = expHkz_22.*phik_2;
    phit_1_store(:,2*i+1) = phik_1;
    phit_2_store(:,2*i+1) = phik_2;
end

%% calculate observable

figure;
set(gcf, 'position', [250 70 1400 900]);
% subplot(2,1,1)
plot(T,sz1sz2_t);

toc;