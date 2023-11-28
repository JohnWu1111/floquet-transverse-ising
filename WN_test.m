clear;
clc;
format long
tic;

%% paramter
L = 1e4;
J = 1;
g = 1;
% k = (2/L:2/L:1)';
k = (1/L:2/L:(L-1)/L)';
dk = 2*pi/L;
x = (0:L)';
xx = (1:L)';
m = 2;
n = 2;

yk = -sinpi(n*k) - J*sinpi(m*k);
zk = -cospi(n*k) + J*cospi(m*k) - g;
dy = (circshift(yk,-1) - circshift(yk,1))/2;
dz = (circshift(zk,-1) - circshift(zk,1))/2;
r2 = yk.^2+zk.^2;
temp = 2*(zk.*dy - yk.*dz)./r2/(2*pi);
WN = sum(temp);

toc;