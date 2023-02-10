%% Setting paths
addpath('module');
alpha = 1;
beta = 1;
param.ne = 300;
param.ni = 100;
param.M        = 100*beta;
param.Mr       = 66*beta;
p = 0.8;
param.p_ee     = p;
param.p_ie     = p;
param.p_ei     = p;
param.p_ii     = p;
param.s_ee     = 5*0.15/alpha*beta/p;
param.s_ie     = 2*0.5/alpha*beta/p;
param.s_ei     = 4.91*0.424/alpha*beta/p;
param.s_ii     = 4.91*0.40/alpha*beta/p;

param.lambda_e = 7000*beta;
param.lambda_i = 7000*beta;
param.tau_ee   = 1.4;
param.tau_ie   = 1.2;
param.tau_ei    = 4.5;
param.tau_ii    = 4.5;
param.tau_re=0;
param.tau_ri=0;
param.duration = 3;
param.delta_time = 0.1;

param.factor_Leak=0;
param.LeakE = 20;
param.s_exe = 1;
param.s_exi = 1;
param.LeakI = 16.7;
param.factor_Leak = inf;
param.ns_ee=alpha;
param.ns_ie=alpha;
param.ns_ei=alpha;
param.ns_ii=alpha;
param.gridsize=0.1;

%%
beta = 1;

param1 = param;
param1.p_ee     = 0.15;
param1.p_ie     = 0.5;
param1.p_ei     = 0.42;
param1.p_ii     = 0.4;
param1.s_ee     = 5;
param1.s_ie     = 2;
param1.s_ei     = 4.91;
param1.s_ii     = 4.91;

tic;
res_lif=model_LIF_FC(param1,[]);
toc;

%%
param2 = param;
param2.s_ei =  4.91*0.414;
tic;
res_lif2=model_LIF_FC(param2,[]);
toc;

%%
figure;
subplot(2,1,1);
rasterplot(res_lif, param);
hold on;
xline(res_lif.MFE_time(1:res_lif.wave_count,1),'g','LineWidth',1);
hold on;
xline(res_lif.MFE_time(1:res_lif.wave_count,2),'--g','LineWidth',1);
xlim([2860, 2920]);
%xlabel('Time (ms)')
set(gca,'fontsize',12,'fontname','Arial');

subplot(2,1,2);
rasterplot(res_lif2, param);
hold on;
xline(res_lif2.MFE_time(1:res_lif4.wave_count,1),'g','LineWidth',1);
hold on;
xline(res_lif2.MFE_time(1:res_lif4.wave_count,2),'--g','LineWidth',1);
xlim([2820, 2900]);
%xlabel('Time (ms)')
set(gca,'fontsize',12,'fontname','Arial');
set(gcf,'Position',[10,10,200,600]);
exportgraphics(gcf,'2-a-1.eps','ContentType','vector')

%% Panel C 

t = 2583;
res = res_lif2;
VE = res_lif2.VE;
VI = res_lif2.VI;
VE_t = VE(t,:)/100;
VI_t = VI(t,:)/100;
edges= -66:2:100;
edges= edges/100;

% Plots
figure;
subplot(2,1,1);
histogram(VE_t,edges,'Normalization', 'probability','FaceColor','red');
hold on;
xline(mean(VE_t),'Color','Green','LineWidth',1);
hold on;
xline(mean(VI_t),'Color','Green','LineWidth',1);
set(gca,'fontsize',15,'fontname','Arial');
ylabel("Excitatory");
subplot(2,1,2);
histogram(VI_t,edges,'Normalization', 'probability','FaceColor','blue');
hold on;
xline(mean(VI_t),'Color','Green','LineWidth',1);
set(gca,'fontsize',15,'fontname','Arial');
ylabel("Inhibitory");
set(gcf,'Position',[10,10,1000,300]);