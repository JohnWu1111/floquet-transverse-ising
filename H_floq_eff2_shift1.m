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
dt = 10;
t0 = 0.5*dt;

Hop_space = zeros(L+1,1);
Gen_space = zeros(L+1,1);
Eli_space = zeros(L+1,1);

%% constructing evolution operator
H1_store = cell(L/2,1);
H2_store = cell(L/2,1);
H_eff_store = cell(L/2,1);

ck = cospi(k);
sk = sinpi(k);
for i = 1:L/2 
    H1 = 2*[g+J*ck(i) -1i*J*sk(i);
                1i*J*sk(i) -(g+J*ck(i))];
    H2 = 2*[-g+J*ck(i) -1i*J*sk(i);
                1i*J*sk(i) -(-g+J*ck(i))];

    expH1 = expm(-1i*t0*H1);
    expH2 = expm(-1i*dt*H2);
    expH3 = expm(-1i*(dt-t0)*H1);
    expH = expH3*expH2*expH1;

    H_eff = 1i*logm(expH)/(2*dt);

    Hop_space = Hop_space + 2*(cospi(k(i)*x)*H_eff(1,1))/L;
    Gen_space = Gen_space - 2i*sinpi(k(i)*x)*H_eff(1,2)/L;
    Eli_space = Eli_space - 2i*sinpi(k(i)*x)*H_eff(2,1)/L;
    
    H1_store{i} = H1;
    H2_store{i} = H2;
    H_eff_store{i} = H_eff;
end

order = sum(abs(Hop_space(1:L/2).*(1:L/2)'))/sum(abs(Hop_space(1:L/2)));

figure
set(gcf, 'position', [250 70 1400 900]);
subplot(2,1,1)
plot(x,abs(Hop_space))
subplot(2,1,2)
plot(x,abs(Gen_space))

toc;