clear;
clc;
format long
tic;

%% paramter
myseed = 1;
rng(myseed)

L = 10000;
dt = 1;
Tmax = 100;
T = 0:dt:Tmax;
nT = length(T);
len = 2^L;
g0 = 0;
J = -1;
g = 0.0;
dk = 2*pi/L;
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

phik_1 = phik0_1;
phik_2 = phik0_2;
phit_1_store = zeros(L/2,nT);
phit_2_store = zeros(L/2,nT);
phit_1_store(:,1) = phik0_1;
phit_2_store(:,1) = phik0_2;

% sz1sz2 = zeros(1,nT);
% Cr = zeros(L,nT);
% Fr = zeros(L,nT);
% n = zeros(1,nT);
% n1n2 = zeros(1,nT);
Gkk = zeros(L,nT);
% n(1) = 2*(phik_1'*phik_1)/L;
% Cr(:,1) = 2*(abs(phik_1).^2)'*cospi(k*x')/L;
% Fr(:,1) = -2i*(phik_1.*conj(phik_2))'*sinpi(k*x')/L;
% n1n2(1) = n(1)^2 + Fr(2)^2 - Cr(2)^2;
% sz1sz2(1) = 1-4*n(1)+4*n1n2(1);
phik_1_all = [phik_1;flip(phik_1)];
phik_2_all = [phik_2;flip(phik_2)];
dphik_1 = (circshift(phik_1_all,-1) - circshift(phik_1_all,1))/2;
dphik_2 = (circshift(phik_2_all,-1) - circshift(phik_2_all,1))/2;
Gkk(:,1) = sqrt(abs(dphik_1).^2+abs(dphik_2).^2-abs(conj(dphik_1).*phik_1_all+conj(dphik_2).*phik_2_all).^2);

for i = 2:nT
    phik_1t = expHkp_11.*phik_1 + expHkp_12.*phik_2;
    phik_2 = expHkp_21.*phik_1 + expHkp_22.*phik_2;
    phik_1 = phik_1t;
    phit_1_store(:,i) = phik_1;
    phit_2_store(:,i) = phik_2;

%     n(i) = 2*(phik_1'*phik_1)/L;
%     Cr(:,i) = 2*(abs(phik_1).^2)'*cospi(k*x')/L;
%     Fr(:,i) = -2i*(phik_1.*conj(phik_2))'*sinpi(k*x')/L;
%     n1n2(i) = n(i)^2 + Fr(2,i)*conj(Fr(2,i)) - Cr(2,i)^2;
%     sz1sz2(i) = 1-4*n(i)+4*n1n2(i);
    phik_1_all = [phik_1;flip(phik_1)];
    phik_2_all = [phik_2;flip(phik_2)];
    dphik_1 = (circshift(phik_1_all,-1) - circshift(phik_1_all,1))/2;
    dphik_2 = (circshift(phik_2_all,-1) - circshift(phik_2_all,1))/2;
    Gkk(:,i) = sqrt(abs(dphik_1).^2+abs(dphik_2).^2-abs(conj(dphik_1).*phik_1_all+conj(dphik_2).*phik_2_all).^2);
end

QV = sum(Gkk);

% collect = [real(phik_1_all);imag(phik_1_all);real(phik_2_all);imag(phik_2_all)];
% figure
% histogram(collect,'NumBins',1000)

% figure
% subplot(1,2,1)
% plot(k,real(phit_1_store(:,end)),k,imag(phit_1_store(:,end)))
% subplot(1,2,2)
% plot(k,real(phit_2_store(:,end)),k,imag(phit_2_store(:,end)))

%% calculate observable

figure;
% set(gcf, 'position', [250 70 1400 900]);
% % subplot(2,1,1)
plot(T,QV);

toc;
