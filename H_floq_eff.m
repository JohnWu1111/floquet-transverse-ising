% calculation of effective floquet Hamiltonian of KICKED Ising model

clear;
clc;
format long
tic;

%% paramter

L = 1000;
J = 0.3; % *pi
g = 0.3; % *pi
k = (2/L:2/L:1)';
x = (0:L)';
xx = (1:L)';

fact = exp(-2i*g);
cc = cospi(2*J);
ss = sinpi(2*J);
expHk_store = cell(L/2,1);
Hk_store = cell(L/2,1);
Hop_space = zeros(L+1,1);
Gen_space = zeros(L+1,1);
Eli_space = zeros(L+1,1);
watch = zeros(L/2,1);
for i = 1:L/2
    expHk = zeros(2);
    expHk(1,1) = fact*(cc + 1i*ss*cospi(k(i)));
    expHk(2,2) = conj(expHk(1,1));
%     expHk(2,2) = 1/fact*(cc - 1i*ss*cospi(k(i)));
    expHk(1,2) = -fact*ss*sinpi(k(i));
    expHk(2,1) = -conj(expHk(1,2));
%     expHk(1,2) = 1/fact*ss*sinpi(k(i));
    expHk_store{i} = expHk;
    Hk = -1i*logm(expHk);
    watch(i) = imag(Hk(1,2));
    Hk_store{i} = Hk;
%     Hop_space = Hop_space + (exp(1i*k(i)*x)*Hk(1,1)-exp(-1i*k(i)*x)*Hk(2,2))/L;
    Hop_space = Hop_space + 2*(cospi(k(i)*x)*Hk(1,1))/L;
    Gen_space = Gen_space - 2i*sinpi(k(i)*x)*Hk(1,2)/L;
    Eli_space = Eli_space - 2i*sinpi(k(i)*x)*Hk(2,1)/L;
end

figure
set(gcf, 'position', [250 70 1400 900]);
subplot(2,1,1)
plot(x,real(Hop_space))
subplot(2,1,2)
plot(x,real(Gen_space))
% subplot(2,1,1)
% plot(x-L/2,circshift(real(Hop_space),L/2))
% subplot(2,1,2)
% plot(x-L/2,circshift(real(Gen_space),L/2))

toc;