% calculation of effective floquet Hamiltonian of transverse Ising model
% with alternating transverse field

clear;
clc;
format long
tic;

%% paramter
L = 1e6;
J = -1;
g = 0;
% k = (2/L:2/L:1)';
k = (1/L:2/L:(L-1)/L+1)';
dk = 2*pi/L;
x = (0:L)';
xx = (1:L)';
% dt_all = 0.001:0.001:1;
% dt0 = 1.1181;
dt0 = 1;
% dt = dt0*pi/(2*sqrt(J^2+g1^2));
dt = dt0*pi/(2*J);
% dt = 100;

%% constructing evolution operator
ck = cospi(k);
sk = sinpi(k);

H_11 = 2*(g+J*ck);
H_12 = 2i*J*sk;

spe = zeros(L,2);
spe(:,1) = -sqrt(H_11.^2+abs(H_12).^2);
spe(:,2) = -spe(:,1);

yk = imag(H_12);
zk = H_11;

dy = (circshift(yk,-1) - circshift(yk,1))/2;
dz = (circshift(zk,-1) - circshift(zk,1))/2;
r2 = yk.^2 + zk.^2;
temp = (zk.*dy - yk.*dz)./r2/(2*pi);
WN = sum(temp);

ftitle = strcat('L = ', num2str(L),', g = ', num2str(g),', dt0 = ', num2str(dt0));
figure('Name',ftitle);
set(gcf, 'position', [250 70 1400 900]);
plot(yk,zk)
% plot3([yk;-yk],[zk;-zk],[k;k+1])
xlabel('yk')
ylabel('zk')

% figure('Name',ftitle);
% set(gcf, 'position', [250 70 1400 900]);
% plot(k,spe)
% xlabel('k')
% ylabel('E_k')


toc;