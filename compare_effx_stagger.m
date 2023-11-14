% calculation of effective floquet Hamiltonian of transverse Ising model
% with alternating transverse field

clear;
clc;
format long
tic;

%% paramter

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
x = (0:L)';
xx = (1:L)';

%% constructing evolution operator
ck = cospi(k);
sk = sinpi(k);



fact = sqrt(g^2+J^2+2*g*J*ck);
sf = sin(2*fact*dt)./fact;

expHkp_11 = cos(2*dt*fact) - 1i*(g+J*ck).*sf;
expHkp_22 = conj(expHkp_11);
expHkp_12 = -J*sk.*sf;
expHkp_21 = -expHkp_12;

fact = sqrt(g^2+J^2-2*g*J*ck);
sf = sin(2*fact*dt)./fact;

expHkm_11 = cos(2*dt*fact) - 1i*(-g+J*ck).*sf;
expHkm_22 = conj(expHkm_11);
expHkm_12 = -J*sk.*sf;
expHkm_21 = -expHkm_12;

expH_11 = expHkm_11.*expHkp_11 + expHkm_12.*expHkp_21;
expH_12 = expHkm_11.*expHkp_12 + expHkm_12.*expHkp_22;
expH_21 = expHkm_21.*expHkp_11 + expHkm_22.*expHkp_21;
expH_22 = expHkm_21.*expHkp_12 + expHkm_22.*expHkp_22;

a = real(expH_11);
b = imag(expH_11);
c = real(expH_12);
d = imag(expH_12);
fact = sqrt(b.^2 + c.^2 + d.^2);

H_eff_11 = 1i*((fact-b).*log(a-1i*fact) + (fact+b).*log(a+1i*fact))./(2*fact);
% H_eff_22 = ((fact+b).*log(a-1i*fact) + (fact-b).*log(a+1i*fact))./(2*fact);
H_eff_22 = -H_eff_11;
H_eff_12 = -(c+d*1i).*(log(a-1i*fact)-log(a+1i*fact))./(2*fact);
H_eff_21 = conj(H_eff_12);

Hop_space = 2*(cospi(x*k')*H_eff_11)/L;
Gen_space = -2i*(sinpi(x*k')*H_eff_12)/L;



% figure
% set(gcf, 'position', [250 70 1400 900]);


toc;