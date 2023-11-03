% calculation of effective floquet Hamiltonian of transverse Ising model
% with alternating transverse field

clear;
clc;
format long
tic;

%% paramter

L = 1000;
J = 1;
g = 1;
% k = (2/L:2/L:1)';
k = (1/L:2/L:(L-1)/L)';
x = (0:L)';
xx = (1:L)';
T = 100;

Hop_space = zeros(L+1,1);
Gen_space = zeros(L+1,1);
Eli_space = zeros(L+1,1);

%% constructing evolution operator
H1_store = cell(L/2,1);
H2_store = cell(L/2,1);
H_eff_store = cell(L/2,1);
for i = 1:L/2
    H1 = 2*[g+J*cospi(k(i)) -1i*J*sinpi(k(i));
                1i*J*sinpi(k(i)) -(g+J*cospi(k(i)))];
    H2 = 2*[-g+J*cospi(k(i)) -1i*J*sinpi(k(i));
                1i*J*sinpi(k(i)) -(-g+J*cospi(k(i)))];

    expH1 = expm(-1i*T*H1);
    expH2 = expm(-1i*T*H2);
    expH = expH2*expH1;

    H_eff = 1i*logm(expH);

    Hop_space = Hop_space + 2*(cospi(k(i)*x)*H_eff(1,1))/L;
    Gen_space = Gen_space - 2i*sinpi(k(i)*x)*H_eff(1,2)/L;
    Eli_space = Eli_space - 2i*sinpi(k(i)*x)*H_eff(2,1)/L;
    
    H1_store{i} = H1;
    H2_store{i} = H2;
    H_eff_store{i} = H_eff;
end

figure
set(gcf, 'position', [250 70 1400 900]);
subplot(2,1,1)
plot(x,real(Hop_space))
subplot(2,1,2)
plot(x,abs(Gen_space))

toc;