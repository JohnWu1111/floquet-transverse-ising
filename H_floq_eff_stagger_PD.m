% calculation of effective floquet Hamiltonian of transverse Ising model
% with alternating transverse field

clear;
clc;
format long
tic;

%% paramter
L_all = 200:200:5000;
% L_all = 1000:1000:5000;
nL = length(L_all);
J = 1;
g = 1;
T_all = 5:5:50;
% T_all = 5:5:100;
nT = length(T_all);

order = zeros(nL,nT);
peak = zeros(nL,nT);

%% constructing evolution operator

for m = 1:nT
    T = T_all(m);
    for n = 1:nL
        L = L_all(n);
        x = (0:L)';
        k = (1/L:2/L:(L-1)/L)';
        ck = cospi(k);
        sk = sinpi(k);
        fact = sqrt(g^2+J^2+2*g*J*ck);
        sf = sin(2*fact*T)./fact;

        expHkp_11 = cos(2*T*fact) - 1i*(g+J*ck).*sf;
        expHkp_22 = conj(expHkp_11);
        expHkp_12 = -J*sk.*sf;
        expHkp_21 = -expHkp_12;

        fact = sqrt(g^2+J^2-2*g*J*ck);
        sf = sin(2*fact*T)./fact;

        expHkm_11 = cos(2*T*fact) - 1i*(-g+J*ck).*sf;
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
        H_eff_12 = -expH_12.*(log(a-1i*fact)-log(a+1i*fact))./(2*fact);
        H_eff_21 = conj(H_eff_12);

%         Hop_space = 2*(cospi(x*k')*H_eff_11)/L;
%         Gen_space = -2i*(sinpi(x*k')*H_eff_12)/L;
        Hop_space = 2*(cos(x*k'*pi)*H_eff_11)/L;
        Gen_space = -2i*(sin(x*k'*pi)*H_eff_12)/L;
        order(n,m) = sum(abs(Hop_space(1:L/2).*(1:L/2)'))/sum(abs(Hop_space(1:L/2)));
        [~,peak(n,m)] = max(abs(Hop_space(1:L/2)));
    end
end

le = cell(1, nT);
for i = 1:nT
    le{i} = strcat('dt = ', num2str(T_all(i)));
end
figure
set(gcf, 'position', [250 70 1400 900]);
plot(L_all,order)
xlabel('L')
ylabel('LR factor')
legend(le)

% le = cell(1, nL);
% for i = 1:nL
%     le{i} = strcat('L = ', num2str(L_all(i)));
% end
% figure
% set(gcf, 'position', [250 70 1400 900]);
% plot(T_all,order)
% xlabel('T')
% ylabel('LR factor')
% legend(le)

toc;