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
g2 = 0.5;
% k = (2/L:2/L:1)';
k = (1/L:2/L:(L-1)/L+1)';
dk = 2*pi/L;
x = (0:L)';
xx = (1:L)';
dt_all = 0.01:0.01:10;
% dt_all = 2;
ndt = length(dt_all);
% dt = 100;

WN1 = zeros(ndt,1);
WN2 = zeros(ndt,1);

%% constructing evolution operator
for n = 1:ndt
    dt = dt_all(n)*pi/(2*sqrt(2));

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

%     H_eff_11 = real(1i*((fact-b).*log(a-1i*fact) + (fact+b).*log(a+1i*fact))./(2*fact));
%     % H_eff_22 = ((fact+b).*log(a-1i*fact) + (fact-b).*log(a+1i*fact))./(2*fact)/(2*T);
%     H_eff_22 = -H_eff_11;
%     H_eff_12 = -(c+d*1i).*(log(a-1i*fact)-log(a+1i*fact))./(2*fact);
%     H_eff_21 = conj(H_eff_12);

    bfact = b./(2*fact);
    lapf = log(a+1i*fact);
    lamf = log(a-1i*fact);
    H_eff_11 = real(1i*((1-bfact).*lamf + (1+bfact).*lapf));
    H_eff_22 = -H_eff_11;
    H_eff_12 = -(c+d*1i).*(lamf-lapf)./(2*fact);
    H_eff_21 = conj(H_eff_12);

%     bfact = b./(2*fact);
%     lp = log(a.^2+fact.^2);
%     lm = log((a+1i*fact)./(a-1i*fact));
%     H_eff_11 = real(1i*(lp+bfact.*lm));
%     H_eff_11 = mod(H_eff_11+pi,2*pi)-pi;
%     H_eff_22 = -H_eff_11;
%     H_eff_12 = (c+d*1i).*lm./(2*fact);
%     H_eff_21 = conj(H_eff_12);

%     spe = zeros(2,L/2);

%     ykr = real(H_eff_12);
%     yki = imag(H_eff_12);
%     ykr = abs(H_eff_12);
%     yki = 0;
    yk = H_eff_12;
    zk = H_eff_11;

%     phi1 = atan((zk+yki)./ykr);
%     phi2 = atan((zk-yki)./ykr);
%     dphi1 = (circshift(phi1,-1) - circshift(phi1,1))/2;
%     dphi2 = (circshift(phi2,-1) - circshift(phi2,1))/2;
%     WN1(n) = sum(dphi1)/(2*pi);
%     WN2(n) = sum(dphi2)/(2*pi);

    phi1 = atan(yk./zk);
    dphi1 = (circshift(phi1,-1) - circshift(phi1,1))/2;
    WN1(n) = sum(dphi1)/(2*pi);
end

ftitle = strcat('L = ', num2str(L),', g1 = ', num2str(g1),', g2 = ', num2str(g2));
figure('Name',ftitle);
set(gcf, 'position', [250 70 1400 900]);
subplot(2,1,1)
plot(dt_all,WN1)
xlabel('/(2*pi/sqrt(2))')
ylabel('Winding Number')
subplot(2,1,2)
plot(dt_all,WN2)
xlabel('/(2*pi/sqrt(2))')
ylabel('Winding Number')

toc;