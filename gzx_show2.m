clear;
clc;
format long
tic;

%% paramter

h0 = 0.0;
h = 0:0.01:2;
nh = length(h0);
dk = 2/1000;
k = dk:dk:1;
nk = length(k);

Gkk_di = 2*(sinpi(k')).^2./(1+h.^2+2*cospi(k')*h);
fact = sqrt(1+h0^2+2*h0*cospi(k'));
temp1 = sinpi(k').^2.*sqrt(((h0+cospi(k')-fact).^2*h.^2.*(h-h0).^2));
temp2 = 1+h.^2+2*cospi(k').*h;
temp3 = cospi(k').*(-2*h0+fact)+h0*(-h0+fact)-1;
Gkk = -2*temp1./temp2./temp3;

Gkk0 = abs(h.*(h-h0)/(h0-1)).*Gkk_di;
Gkk_diff = Gkk-Gkk0;

Gkk = abs(pi-k'*pi).^2./(1+h.^2+2*cospi(k')*h);

QV = 2*sum(Gkk).*dk*pi;

% slope_l = (QV(100)-QV(99))/0.01;
% slope_r = (QV(103)-QV(102))/0.01;
% slope_diff = slope_l - slope_r

figure
% mesh(h,k,Gkk)
% xlabel('h')
% ylabel('k')
plot(h,QV)

toc;