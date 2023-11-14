% calculation of effective floquet Hamiltonian of transverse Ising model
% with alternating transverse field

clear;
clc;
format long
tic;

%% paramter
L = 1e5;
J = 1;
g1 = 1;
g2 = -1;
% k = (2/L:2/L:1)';
k = (1/L:2/L:(L-1)/L)';
dk = 2*pi/L;
x = (0:L)';
xx = (1:L)';
% dt_all = 0.001:0.001:1;
dt = 5*pi/(2*sqrt(2));
% dt = 100;

%% constructing evolution operator

ck = cospi(k);
sk = sinpi(k);
fact = sqrt(g1^2+J^2+2*g1*J*ck);
sf = sin(2*fact*dt)./fact;

expHkp_11 = cos(2*dt*fact) - 1i*(g1+J*ck).*sf;
expHkp_22 = conj(expHkp_11);
expHkp_12 = -J*sk.*sf;
expHkp_21 = -expHkp_12;

fact = sqrt(g2^2+J^2+2*g2*J*ck);
sf = sin(2*fact*dt)./fact;

expHkm_11 = cos(2*dt*fact) - 1i*(g2+J*ck).*sf;
expHkm_22 = conj(expHkm_11);
expHkm_12 = -J*sk.*sf;
expHkm_21 = -expHkm_12;

% fact = sqrt(g1^2+J^2+2*g1*J*ck);
% sf = sin(2*fact*(dt-t0))./fact;
%
% expHk2_11 = cos(2*(dt-t0)*fact) - 1i*(g1+J*ck).*sf;
% expHk2_22 = conj(expHk2_11);
% expHk2_12 = -J*sk.*sf;
% expHk2_21 = -expHk2_12;

expH_11 = expHkm_11.*expHkp_11 + expHkm_12.*expHkp_21;
expH_12 = expHkm_11.*expHkp_12 + expHkm_12.*expHkp_22;
expH_21 = expHkm_21.*expHkp_11 + expHkm_22.*expHkp_21;
expH_22 = expHkm_21.*expHkp_12 + expHkm_22.*expHkp_22;

% expH_11 = expHk2_11.*expH_11 + expHk2_12.*expH_21;
% expH_12 = expHk2_11.*expH_12 + expHk2_12.*expH_22;
% expH_21 = expHk2_21.*expH_11 + expHk2_22.*expH_21;
% expH_22 = expHk2_21.*expH_12 + expHk2_22.*expH_22;

a = real(expH_11);
b = imag(expH_11);
c = real(expH_12);
d = imag(expH_12);
fact = sqrt(b.^2 + c.^2 + d.^2);

H_eff_11 = real(1i*((fact-b).*log(a-1i*fact) + (fact+b).*log(a+1i*fact))./(2*fact));
% H_eff_22 = ((fact+b).*log(a-1i*fact) + (fact-b).*log(a+1i*fact))./(2*fact)/(2*T);
H_eff_22 = -H_eff_11;
H_eff_12 = -(c+d*1i).*(log(a-1i*fact)-log(a+1i*fact))./(2*fact);
H_eff_21 = conj(H_eff_12);

spe = zeros(2,L/2);

zk = H_eff_11;
yk = imag(H_eff_12);

% dy = (yk - circshift(yk,-1));
% dz = (zk - circshift(zk,-1));
dy = (circshift(yk,1) - circshift(yk,-1))/2;
dz = (circshift(zk,1) - circshift(zk,-1))/2;
r2 = yk.^2+zk.^2;
temp = 2*(zk.*dy - yk.*dz)./r2/(2*pi);
WN = sum(temp);

ftitle = strcat('L = ', num2str(L),', g1 = ', num2str(g1),', g2 = ', num2str(g2));
figure('Name',ftitle);
set(gcf, 'position', [250 70 1400 900]);
plot([yk;-yk],[zk;-zk])
xlabel('yk')
ylabel('zk')


toc;