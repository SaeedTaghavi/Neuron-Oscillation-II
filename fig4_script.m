  %% Setting Path
addpath('module');
%% Parameters for Model LIF

alpha = 1;
beta = 1;
param.ne = 300;
param.ni = 100;
param.M        = 100*beta;
param.Mr       = 66*beta;
% param.p_ee     = 0.15;
% param.p_ie     = 0.5;
% param.p_ei     = 0.415;
% param.p_ii     = 0.4;
param.p_ee     = 0.8;
param.p_ie     = 0.8;
param.p_ei     = 0.8;
param.p_ii     = 0.8;
% param.s_ee     = 5/alpha*beta;
% param.s_ie     = 2/alpha*beta;
% param.s_ei     = 4.91/alpha*beta;
% param.s_ii     = 4.91/alpha*beta;
param.s_ee     = 5*0.15/alpha*beta/param.p_ee;
param.s_ie     = 2*0.5/alpha*beta/param.p_ie;
param.s_ei     = 4.91*0.44  /alpha*beta/param.p_ei;
param.s_ii     = 4.91*0.4/alpha*beta/param.p_ii;
param.s_exe = 1;
param.s_exi = 1;

param.lambda_e = 7000*beta;
param.lambda_i = 7000*beta;
param.tau_ee   = 1.4;
param.tau_ie   = 1.2;
param.tau_ei    = 4.5;
param.tau_ii    = 4.5;
param.tau_re=0;
param.tau_ri=0;
param.duration = 3;

param.factor_Leak=0;
param.LeakE = 20;

param.LeakI = 16.7;

param.factor_Leak = inf;
param.ns_ee=alpha;
param.ns_ie=alpha;
param.ns_ei=alpha;
param.ns_ii=alpha;
param.gridsize=0.1;

% Model LIF
tic;
res_lif_0_001=model_LIF3_fixconnection(param,[]);
toc;
figure;
rasterplot(res_lif_0_001, param);
xlim([2000, 3000]);

%% MFE figure

figure;
rasterplot(res_lif, param);
hold on;
xline(res_lif.MFE_time(1:res_lif.wave_count,1),'b','LineWidth',1);
hold on;
xline(res_lif.MFE_time(1:res_lif.wave_count,2),'LineWidth',1);
xlim([2000, 3000]);

%% 代码记录

%sexe
param2.s_exe = 0.2+(i-1)*0.4;
param2.s_exi = 0.2+(i-1)*0.4;
param2.lambda_e = 7000*beta/param2.s_exe;
param2.lambda_i = 7000*beta/param2.s_exi;
%sei
param2.s_ei=(0.40+i*0.0005)*4.91;
%p
param2.p_ee     = 0.25+0.015/2*i;
param2.p_ie     = 0.25+0.015/2*i;
param2.p_ei     = 0.25+0.015/2*i;
param2.p_ii     = 0.25+0.015/2*i;
param2.s_ee     = 5*0.15/alpha*beta/param2.p_ee;
param2.s_ie     = 2*0.5/alpha*beta/param2.p_ie;
param2.s_ei     = 4.91*0.41/alpha*beta/param2.p_ei;
param2.s_ii     = 4.91*0.40/alpha*beta/param2.p_ii;

param2.p_ee     = 1/(1+0.03*i);
param2.p_ie     = 1/(1+0.03*i);
param2.p_ei     = 1/(1+0.03*i);
param2.p_ii     = 1/(1+0.03*i);
param2.s_ee     = 5*0.15/alpha*beta/param2.p_ee;
param2.s_ie     = 2*0.5/alpha*beta/param2.p_ie;
param2.s_ei     = 4.91*0.424/alpha*beta/param2.p_ei;
param2.s_ii     = 4.91*0.40/alpha*beta/param2.p_ii;

%tau
param2.tau_re=i*0.2;
param2.tau_ri=i*0.2;

%tau_e
param2.tau_ee   = 1.4*(0.5+i*0.05);
param2.tau_ie   = 1.2*(0.5+i*0.05);
%%
dmean_v_exnoise=zeros(50,600);

parfor i=1:50
    param2=param;
    i

param2.tau_re=i*0.1;
param2.tau_ri=i*0.1;
    res_lif=model_LIF3(param2,[]);
    ind=find(res_lif.MFE_time(:,1)~=0,1,"last");
    MFE_time=res_lif.MFE_time(1:ind,1);
    dm=mean(res_lif.VE(ceil(MFE_time*10-20),:),2)-mean(res_lif.VI(ceil(MFE_time*10-20),:),2);
    dm=clean_dm(dm);
    dm=dm(ceil(length(dm)*0.3):end)';
    dmean_v_exnoise(i,:)=[length(dm)-1,dm(2:end)-dm(1:end-1),zeros(1,600-length(dm))];
end
%%
dens=gauss_density(dmean_v_exnoise,[40:-1:-40]);
imagesc(dens);
colormap turbo;
caxis([0 0.12]);
xticks([1 10 20 30 40 50]);

xticklabels([0 1 2 3 4 5]);
yticks([1 21 41 61 81]);
yticklabels([0.4 0.2 0 -0.2 -0.4])
xlabel('');
ylabel('\Deltam');
set(gcf,'Position',[10,10,1200,600]);
set(gca,'FontSize',30);
set(gca,'FontName','Arial');
%%
video(res_lif_08)

%%
function dm=clean_dm(dm)
l=length(dm);
i=1;
while i<l
    dm(i)=(dm(i)+dm(i+1)<50)*dm(i);
    i=i+1;
end
dm=dm(dm~=0);
end


%%
function []=video(res_lif)
for i=300:3:1000
subplot(2,1,1)
histogram(res_lif.VE(i,:),[-66:1:100]);
ylim([0 30]);
subplot(2,1,2)
histogram(res_lif.VI(i,:),[-66:1:100]);
ylim([0 30]);
pause(0.6);
end
end