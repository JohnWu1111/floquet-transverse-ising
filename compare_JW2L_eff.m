clear;
clc;
format long
tic;

%% paramter

L = 4;
num_T = 200;
T = 0:num_T;
nT = num_T+1;
dt = 1;
len = 2^L;
g0 = 0.1;
J = -1;
g = 1;
% k = (2/L:2/L:1)';
k = (1/L:2/L:(L-1)/L)';
x = (0:L)';
xx = (0:L-1)';

%% computation of H_eff
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

H_eff_11 = 1i*((fact-b).*log(a-1i*fact) + (fact+b).*log(a+1i*fact))./(2*fact)/(2*dt);
% H_eff_22 = ((fact+b).*log(a-1i*fact) + (fact-b).*log(a+1i*fact))./(2*fact);
H_eff_22 = -H_eff_11;
H_eff_12 = -expH_12.*(log(a-1i*fact)-log(a+1i*fact))./(2*fact)/(2*dt);
H_eff_21 = conj(H_eff_12);

Hop_space = real(2*(cospi(xx*k')*H_eff_11))/L;
Gen_space = 2i*(sinpi(xx*k')*H_eff_12)/L;

order = sum(abs(Hop_space(1:L/2).*(1:L/2)'))/sum(abs(Hop_space(1:L/2)));

%% construction of Hamiltonian
A = zeros(L);
for i = 1:L-1
    A(i,i+1) = 1/2;
    A(i+1,i) = 1/2;
end
A(L,1) = -1/2;
A(1,L) = -1/2;

B = zeros(L);
for i = 1:L-1
    B(i,i+1) = 1/2;
    B(i+1,i) = -1/2;
end
B(L,1) = -1/2;
B(1,L) = 1/2;

%% initial state
Hxx = [A,B;
      -B,-A];
Hz = [ones(L,1);-ones(L,1)];
H0 = g0*diag(Hz) + J*Hxx;
[V, D] = eig(H0);
e = diag(D);
e_GS = sum(e(1:L));

V_GS = V(:,1:L);
G0 = V_GS*V_GS';
G = G0;
Gt_store = cell(nT,1);
Gt_store{1} = G;

%% time evolution
sz1sz2 = zeros(1,nT);
sz1sz2(1) = 1-4*G(1,1)+4*(G(1,1)^2+G(1+L,2)^2-G(1,2)^2);
n1 = zeros(1,nT);
n1(1) = G(1,1);

% exp(Hm)*exp(Hp)
[Vp,Dp] = eig(J*Hxx+g*diag(Hz));
ep = diag(Dp);
[Vm,Dm] = eig(J*Hxx-g*diag(Hz));
em = diag(Dm);

trans_p = exp(2i*ep);
trans0_p = Vp.*trans_p'*Vp';
trans_m = exp(2i*em);
trans0_m = Vm.*trans_m'*Vm';

for i = 2:nT
    % Hp
    G = trans0_p'*G*trans0_p;

    % Hm
    G = trans0_m'*G*trans0_m;
    Gt_store{i} = G;
    sz1sz2(i) = 1-4*G(1,1)+4*(G(1,1)^2+abs(G(1+L,2))^2-G(1,2)^2);
    n1(i) = G(1,1);
end

% % exp(-1i*H_eff)
% H_hop = zeros(L);
% H_gen = zeros(L);
% temp_Hop = Hop_space;
% temp_Gen = Gen_space;
% for i = 1:L
%     H_hop(:,i) = temp_Hop;
%     H_gen(:,i) = temp_Gen;
%     temp_Hop = [-temp_Hop(L);temp_Hop(1:L-1)];
%     temp_Gen = [-temp_Gen(L);temp_Gen(1:L-1)];
% end
% 
% % H_eff = [H_hop,H_gen;H_gen',-H_hop];
% H_eff = [H_hop/2,H_gen/2;-conj(H_gen)/2,-H_hop/2];
% [V,D] = eig(H_eff);
% e = diag(D);
% 
% trans = exp(2i*e*2*dt);
% trans0 = V.*trans'*V';
% 
% for i = 2:nT
%     G = trans0'*G*trans0;
%     Gt_store{i} = G;
%     sz1sz2(i) = 1-4*G(1,1)+4*(G(1,1)^2+abs(G(1+L,2))^2-G(1,2)^2);
%     n1(i) = G(1,1);
% end

%% calculate observable

figure;
set(gcf, 'position', [250 70 1400 900]);
% subplot(2,1,1)
plot(T,sz1sz2);

toc;