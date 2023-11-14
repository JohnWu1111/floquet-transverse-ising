% calculation of effective floquet Hamiltonian of transverse Ising model
% with alternating transverse field

clear;
clc;
format long
tic;

%% paramter

L = 2000;
J = 1;
g1 = 1;
g2 = 2;
% k = (2/L:2/L:1)';
k = (1/L:2/L:(L-1)/L)';
x = (0:L)';
xx = (1:L)';
T = 10;

%% constructing evolution operator
ck = cospi(k);
sk = sinpi(k);
fact = sqrt(g1^2+J^2+2*g1*J*ck);
sf = sin(2*fact*T)./fact;

expHkp_11 = cos(2*T*fact) - 1i*(g1+J*ck).*sf;
expHkp_22 = conj(expHkp_11);
expHkp_12 = -J*sk.*sf;
expHkp_21 = -expHkp_12;

fact = sqrt(g2^2+J^2+2*g2*J*ck);
sf = sin(2*fact*T)./fact;

expHkm_11 = cos(2*T*fact) - 1i*(g2+J*ck).*sf;
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

H_eff_11 = 1i*((fact-b).*log(a-1i*fact) + (fact+b).*log(a+1i*fact))./(2*fact)/(2*T);
% H_eff_22 = ((fact+b).*log(a-1i*fact) + (fact-b).*log(a+1i*fact))./(2*fact)/(2*T);
H_eff_22 = -H_eff_11;
H_eff_12 = -(c+d*1i).*(log(a-1i*fact)-log(a+1i*fact))./(2*fact)/(2*T);
H_eff_21 = conj(H_eff_12);

Hop_space = 2*(cospi(x*k')*H_eff_11)/L;
Gen_space = -2i*(sinpi(x*k')*H_eff_12)/L;

order = sum(abs(Hop_space(1:L/2).*(1:L/2)'))/sum(abs(Hop_space(1:L/2)));
[~,peak] = max(abs(Hop_space(1:L/2)));

figure
set(gcf, 'position', [250 70 1400 900]);
subplot(2,1,1)
plot(x,abs(Hop_space))
subplot(2,1,2)
plot(x,abs(Gen_space))

fit_x = 1:L/2;
fit_y = real(Gen_space(2:L/2+1))';

toc;