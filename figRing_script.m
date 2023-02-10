%% Setting Path

addpath('module');

%% Parameters


alpha = 1;
beta = 1;
param.ne = 300;
param.ni = 100;
param.M        = 100*beta;
param.Mr       = 66*beta;
param.tau_ee   = 1.4;
param.tau_ie   = 1.2;
param.tau_ei    = 4.5;
param.tau_ii    = 4.5;
param.tau_re=0;
param.tau_ri=0;
param.duration = 3;
param.delta_time = 0.1;
param.s_exe = 1;
param.s_exi = 1;
param.lambda_e = 7000*beta;
param.lambda_i = 7000*beta;
param.gridsize=0.1;


%% Models 


p = 0.2;
param2 = param;
param2.p = p;
param2.p_ee     = p;
param2.p_ie     = p;
param2.p_ei     = p;
param2.p_ii     = p;
param2.s_ee     = 5*0.15/alpha*beta/p;
param2.s_ie     = 2*0.5/alpha*beta/p;
param2.s_ei     = 4.91*0.421/alpha*beta/p;
param2.s_ii     = 4.91*0.40/alpha*beta/p;


% Model LIF

tic;
res_lif=model_LIF_FC(param2,[]);
toc;

% Model Ring

tic;
res_ring = model_LIF_Ring(param2,[]);
toc;

%% 
figure;
subplot(2,1,1);
rasterplot(res_lif, param2);
box on;
xlabel('Time (ms)');
xlim([2000,3000]);
set(gca,'fontsize',15,'fontname','Arial');
subplot(2,1,2);
rasterplot(res_ring, param2);
box on;
xlabel('Time (ms)');
xlim([2000,3000]);
set(gca,'fontsize',15,'fontname','Arial');
set(gcf,'Position',[10,10,1200,600]);