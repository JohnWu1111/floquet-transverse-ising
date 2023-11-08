clear;
clc;
format long
tic;

%% paramter
myseed = 1;
rng(myseed)

L = 4;
num_T = 200;
T = 0:2*num_T;
nT = 2*num_T+1;
dt = 1;
len = 2^L;
g0 = 0.1;
J = -1;
g = 1;
k = (0:2/L:2-2/L)';
x = (0:L-1)';
xx = (1:L)';

phi0 = 