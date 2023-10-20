% Since \epsilon is small, we can ignore terms o(\epsilon^5) in the
% evolution operator(matrix)

clear;
clc;
format long
tic;

%% paramter
myseed = 1;
rng(myseed)

L = 18;
num_T = 100;
T = 0:2*num_T;
nT = 2*num_T+1;
len = 2^L;
g0 = pi/2+0.1;
J0 = pi/4;
it = floor(L/2);

sigmaz = [1;-1];
sigmax = sparse([0 1;1 0]);
I2 = eye(2);
II2 = ones(2);

%% construction of Hamiltonian and observable
Hz = zeros(len,1);
for i = 1:L-1
    Hz_temp = ones(2^(i-1),1);
    Hz_temp = kron(Hz_temp, sigmaz);
    Hz_temp = kron(Hz_temp, sigmaz);
    Hz_temp = kron(Hz_temp, ones(2^(L-i-1),1));
    Hz = Hz + J0*Hz_temp;
end
% PBC
Hz_temp = sigmaz;
Hz_temp = kron(Hz_temp, ones(2^(L-2),1));
Hz_temp = kron(Hz_temp, sigmaz);
Hz = Hz + J0*Hz_temp;

expHx1 = I2;
expHx2 = eye(4);
expHx3 = eye(8);
% expHx4 = eye(16);
expHx1 = kron(sigmax, expHx1) + kron(I2, sparse(1:2,2:-1:1,1));

expHx2 = kron(sigmax, expHx2) + kron(I2, expHx1);
expHx1 = kron(sigmax, expHx1) + kron(I2, sparse(1:4,4:-1:1,1));

expHx3 = kron(sigmax, expHx3) + kron(I2, expHx2);
expHx2 = kron(sigmax, expHx2) + kron(I2, expHx1);
expHx1 = kron(sigmax, expHx1) + kron(I2, sparse(1:8,8:-1:1,1));

for i = 5:L
%     expHx4 = kron(sigmax, expHx4) + kron(I2, expHx3);
    expHx3 = kron(sigmax, expHx3) + kron(I2, expHx2);
    expHx2 = kron(sigmax, expHx2) + kron(I2, expHx1);
    expHx1 = kron(sigmax, expHx1) + kron(I2, sparse(1:2^(i-1),2^(i-1):-1:1,1));   
end

[row1,col1] = find(expHx1);
[row2,col2] = find(expHx2);
[row3,col3] = find(expHx3);
% [row4,col4] = find(expHx4);
expHx1 = sparse(row1,col1,(-1i*sin(g0))^(L-1)*cos(g0));
expHx2 = sparse(row2,col2,(-1i*sin(g0))^(L-2)*(cos(g0))^2);
expHx3 = sparse(row3,col3,(-1i*sin(g0))^(L-3)*(cos(g0))^3);
% expHx4 = sparse(row4,col4,(-1i*sin(g0))^(L-4)*(cos(g0))^4);
delete row1 row2 row3 col1 col2 col3
expHxL = sparse(1:len,len:-1:1,(-1i*sin(g0))^L);
expHx = expHx1 + expHx2 + expHx3 + expHxL;
% expHx = expHx1 + expHx2 + expHx3 + expHx4 + expHxL;

% temp = 1;
% for i = 1:L
%     temp = temp*(-sin(g(i)));
% end
% expHx = sparse(1:len,len:-1:1,temp);

matrix_sz = zeros(len,L);
for i = 1:L
    sz_temp = ones(2^(i-1),1);
    sz_temp = kron(sz_temp,sigmaz);
    sz_temp = kron(sz_temp,ones(2^(L-i),1));
    matrix_sz(:,i) = sz_temp;
end

%% time evolution

% phi0 = rand(len,1);
% phi0 = phi0./sqrt(sum(phi0.^2));

[~,index0] = min(Hz);
phi0 = zeros(len,1);
pos = index0(1);

sz_GS0 = zeros(L,1);
temp = pos-1;
for i = 1:L
    sz_GS0(L-i+1) = mod(temp,2);
    temp = floor(temp/2);
end
sz_GS = 2*sz_GS0-1;
phi0(pos) = 1;

% flip a single point
sz_GS_new = sz_GS;
sz_GS_new(it) = - sz_GS(it);
sz_GS0_new = (sz_GS_new+1)/2;
temp = 0;
for i = 1:L
    temp = temp*2;
    temp = temp + sz_GS0_new(i);    
end
pos_new = temp+1;
phi0_e = zeros(len,1);
phi0_e(pos_new) = 1;
phi0_e = phi0_e./sqrt(sum(phi0_e.^2));

phi = phi0;
phit_store = zeros(len,nT);
phit_store(:,1) = phi0;
phi_e = phi0_e;
phit_e_store = zeros(len,nT);
phit_e_store(:,1) = phi0_e;

trans_z = exp(-1i*Hz);

now = 2;
for i = 1:num_T
    % Hz
    phi = trans_z.*phi;
    phit_store(:,now) = phi;
    phi_e = trans_z.*phi_e;
    phit_e_store(:,now) = phi_e;
    now = now + 1;

    % Hx
    phi = expHx*phi;
    phit_store(:,now) = phi;

    phi_e = expHx*phi_e;
    phit_e_store(:,now) = phi_e;
    now = now + 1;
end

%% calculate observable

sz = zeros(L,nT);
sz_e = zeros(L,nT);
for i = 1:L
    sz(i,:) = sum(conj(phit_store).*(matrix_sz(:,i).*phit_store));
    sz_e(i,:) = sum(conj(phit_e_store).*(matrix_sz(:,i).*phit_e_store));
end
sz_diff = sz_e - sz;

sz_mean = sum(sz_GS.*sz)/L;
sz_e_mean = sum(sz_GS.*sz_e)/L;

x = (1:L)';
sz_diff_std = sqrt(sum(abs(sz_diff).*(x-it).^2));

figure;
set(gcf, 'position', [250 70 1400 900]);
subplot(1,2,1)
plot(T,sz_diff_std)

subplot(1,2,2)
% plot(log(T(1:50)),log(sz_diff_std(1:50)));
% plot(T,sz_e_mean);
imagesc(sz_diff)

toc;

% function y = kron_next(A)
%     la = length(A);
%     [rowA,colA] = find(A);
% %         y((rowA-1)*2+1,(colA-1)*2+2) = 1;
% %         y((rowA-1)*2+2,(colA-1)*2+1) = 1;
%     y1 = sparse((2-1)*la+rowA, (1-1)*la+colA, 1, la*2, la*2);
%     y2 = sparse((1-1)*la+rowA, (2-1)*la+colA, 1, la*2, la*2);
% 
%     y3 = sparse((1-1)*la+(1:la), (1-1)*la+(la:-1:1), 1, la*2, la*2);
%     y4 = sparse((2-1)*la+(1:la), (2-1)*la+(la:-1:1), 1, la*2, la*2);
% 
%     y = y1 + y2 + y3 + y4;
% end
