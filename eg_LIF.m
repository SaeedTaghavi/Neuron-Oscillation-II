%% Setting Path

addpath('module');
load('param.mat');

%% Parameters for Model LIF

beta = 10;
p = 0.8;
param.p_ee = p;
param.p_ei = p;
param.p_ii = p;
param.p_ie = p;
param.duration =  3;
alpha = 1;

param3 = param;
param3.M        = 100*beta;
param3.Mr       = 66*beta;
% param3.s_ee     = 5*0.15/alpha*beta/p*300/299;
% param3.s_ie     = 2*0.5/alpha*beta/p*300/299;
% param3.s_ei     = 4.91*0.421/alpha*beta/p*100/99;
% param3.s_ii     = 4.91*0.40/alpha*beta/p*100/99;
param3.s_ee     = 5*0.15/alpha*beta/p;
param3.s_ie     = 2*0.5/alpha*beta/p;
param3.s_ei     = 4.91*0.425/alpha*beta/p;
param3.s_ii     = 4.91*0.40/alpha*beta/p;
param3.tau_ie = 1.2*0.5;
param3.tau_ee = 1.4*0.5;
param3.lambda_e = 7000*beta;
param3.lambda_i = 7000*beta;

tic;
res_lif=model_LIF_FC(param3,[]);
toc;

%%
figure;
rasterplot(res_lif, param3);
xlim([2000,3000]);