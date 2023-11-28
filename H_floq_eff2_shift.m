% calculation of effective floquet Hamiltonian of transverse Ising model
% with alternating transverse field

clear;
clc;
format long
tic;

%% paramter

L = 50;
J = -1;
g1 = 2;
g2 = -2;
% k = (2/L:2/L:1)';
k = (1/L:2/L:(L-1)/L)';
x = (0:L)';
xx = (1:L)';
% dt = 5*pi/sqrt(2);
% dt = 5*pi/2;
dt = 20;
t0 = 0.5*dt;

%% constructing evolution operator
ck = cospi(k);
sk = sinpi(k);
fact = sqrt(g1^2+J^2+2*g1*J*ck);
sf = sin(2*fact*t0)./fact;

expHkp_11 = cos(2*t0*fact) - 1i*(g1+J*ck).*sf;
expHkp_22 = conj(expHkp_11);
expHkp_12 = -J*sk.*sf;
expHkp_21 = -expHkp_12;

fact = sqrt(g2^2+J^2+2*g2*J*ck);
sf = sin(2*fact*dt)./fact;

expHkm_11 = cos(2*dt*fact) - 1i*(g2+J*ck).*sf;
expHkm_22 = conj(expHkm_11);
expHkm_12 = -J*sk.*sf;
expHkm_21 = -expHkm_12;

fact = sqrt(g1^2+J^2+2*g1*J*ck);
sf = sin(2*fact*(dt-t0))./fact;

expHk2_11 = cos(2*(dt-t0)*fact) - 1i*(g1+J*ck).*sf;
expHk2_22 = conj(expHk2_11);
expHk2_12 = -J*sk.*sf;
expHk2_21 = -expHk2_12;

expH_11 = expHkm_11.*expHkp_11 + expHkm_12.*expHkp_21;
expH_12 = expHkm_11.*expHkp_12 + expHkm_12.*expHkp_22;
expH_21 = expHkm_21.*expHkp_11 + expHkm_22.*expHkp_21;
expH_22 = expHkm_21.*expHkp_12 + expHkm_22.*expHkp_22;

expH_11 = expHk2_11.*expH_11 + expHk2_12.*expH_21;
expH_12 = expHk2_11.*expH_12 + expHk2_12.*expH_22;
expH_21 = expHk2_21.*expH_11 + expHk2_22.*expH_21;
expH_22 = expHk2_21.*expH_12 + expHk2_22.*expH_22;

a = real(expH_11);
b = imag(expH_11);
c = real(expH_12);
d = imag(expH_12);
fact = sqrt(b.^2 + c.^2 + d.^2);

H_eff_11 = 1i*((fact-b).*log(a-1i*fact) + (fact+b).*log(a+1i*fact))./(2*fact);
% H_eff_22 = ((fact+b).*log(a-1i*fact) + (fact-b).*log(a+1i*fact))./(2*fact)/(2*T);
H_eff_22 = -H_eff_11;
H_eff_12 = -(c+d*1i).*(log(a-1i*fact)-log(a+1i*fact))./(2*fact);
H_eff_21 = conj(H_eff_12);

Hop_space = 2*(cos(x*k'*pi)*H_eff_11)/L;
Gen_space = -2i*(sin(x*k'*pi)*H_eff_12)/L;

order = sum(abs(Hop_space(1:L/2).*(1:L/2)'))/sum(abs(Hop_space(1:L/2)));
[~,peak] = max(abs(Hop_space(1:L/2)));

ftitle = strcat('L = ', num2str(L),', g1 = ', num2str(g1),', g2 = ', num2str(g2),', dt = ', num2str(dt),', t0 = ', num2str(t0));
figure('Name',ftitle);
set(gcf, 'position', [250 70 1400 900]);
subplot(2,2,1)
loglog(1:L/2,abs(Hop_space(1:L/2)))
subplot(2,2,2)
semilogy(1:L/2,abs(Hop_space(1:L/2)))
subplot(2,2,3)
loglog(1:L/2,abs(Gen_space(1:L/2)))
subplot(2,2,4)
semilogy(1:L/2,abs(Gen_space(1:L/2)))

fit_x = peak-1:L/2;
fit_y = abs(Hop_space(peak:L/2+1))';

toc;