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
param.p_ee     = 1;
param.p_ie     = 1;
param.p_ei     = 1;
param.p_ii     = 1;
% param.s_ee     = 5/alpha*beta;
% param.s_ie     = 2/alpha*beta;
% param.s_ei     = 4.91/alpha*beta;
% param.s_ii     = 4.91/alpha*beta;
param.s_ee     = 5*0.15/alpha*beta;
param.s_ie     = 2*0.5/alpha*beta;
param.s_ei     = 4.91*0.421/alpha*beta;
param.s_ii     = 4.91*0.40/alpha*beta;
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
param.delta_time = 0.1;

param.factor_Leak=0;
param.LeakE = 20;

param.LeakI = 16.7;
param.factor_Leak = inf;
param.ns_ee=alpha;
param.ns_ie=alpha;
param.ns_ei=alpha;
param.ns_ii=alpha;
param.gridsize=0.1;

%% 这段不用跑
param.duration=0.3;
init=struct;
init.ve=sqrt(75)*randn(1,300)+5;
init.vi=sqrt(80)*randn(1,100);
m=max(max(init.ve),max(init.vi))-100;
init.ve=init.ve-m-1;
init.vi=init.vi-m-1;
init.he=zeros(2,300);
init.hi=zeros(2,100);
res_lif=model_LIF3(param,init);

figure;
rasterplot(res_lif, param);
hold on;
for i =1:res_lif.wave_count
    xline(res_lif.MFE_time(i,1),'b','LineWidth',1);
    hold on;
    xline(res_lif.MFE_time(i,2),'LineWidth',1);
    hold on;
end
%% fig3
[res,dm,res_lif,index]=dm_plot(param,'test');

%%
[res,dm,index]=invariant_plot(param,'test');
%%
function [res,dm,res_lif,index]=dm_plot(param,name)
res_lif=model_LIF3(param,[]);

dm=zeros(51,24);
res=cell(51,24);
index=zeros(51,24);

xrange=[-20:30];
parfor d=1:numel(xrange)
    param2=param;
    param2.duration=0.08;
    pardm=zeros(1,24);
    parres=cell(1,24);
    parind=zeros(1,24);
    for trial=1:24
        trial

        init=struct;
        init.ve=sqrt(75)*randn(1,300)+xrange(d);
        init.vi=sqrt(80)*randn(1,100);
        m=max(max(init.ve),max(init.vi))-100;
        init.ve=init.ve-m-1;
        init.vi=init.vi-m-1;
        init.he=zeros(2,300);
        init.hi=zeros(2,100);

        res_lif_small=model_LIF_FC(param2,init);
        ind=find(res_lif_small.wave_spike_count>300);
        ind=min(ind);
        if res_lif_small.MFE_time(ind,1)>10
            pos=ceil(res_lif_small.MFE_time(ind,1)*10-20);
            pardm(trial)=mean(res_lif_small.VE(pos,:))-mean(res_lif_small.VI(pos,:));
        else
            pos=ceil(res_lif_small.MFE_time(ind+1,1)*10-20);
            pardm(trial)=mean(res_lif_small.VE(pos,:))-mean(res_lif_small.VI(pos,:));
        end
        parres{trial}=res_lif_small;
        parind(trial)=pos;
    end
    dm(d,:)=pardm;
    res(d,:)=parres;
    index(d,:)=parind;
end

subplot(2,5,[1,2,6,7]);
xc=xrange'*ones(1,24);
scatter(xc(:),dm(:),'.');
grid on;
hold on;
plot([-30:30],[-30:30]);
xlim([-15 30]);
ylim([-15 30]);
xticks([-30 -25 -20 -15 -10 -5 0 5 10 15 20 25 30])
yticks([-30 -25 -20 -15 -10 -5 0 5 10 15 20 25 30])
xlabel('difference of mean before the N-1th MFE');
ylabel('difference of mean before the Nth MFE');
set(gca,'FontSize',11);
subplot(2,5,[3,4,5]);
rasterplot(res_lif,param);
xlim([2500 3000]);
title('rasterplot');
set(gca,'FontSize',11);
subplot(2,5,[8,9,10]);
dmean_wholetime=mean(res_lif.VE(2:end,:),2)-mean(res_lif.VI(2:end,:),2);
plot(dmean_wholetime(1:10:end));
set(gca,'FontSize',11);
xlim([2500 3000]);
grid on;
yticks([-50 -25 0 25 50]);
xlabel('time(ms)');
title('difference of mean');
set(gcf,'Position',[10,10,1200,420]);

% saveas(gcf ,['outputs/dm mapping graph/parameter_searching/',name,'.fig']);
% close(gcf);
end

%%
function [res,dm,index]=invariant_plot(param,name)


dm=zeros(1,300);
res=cell(1,300);
index=zeros(1,300);



for trial=2:length(dm)
    param2=param;
    param2.duration=0.08;
    param2.s_ei     = 4.91*0.41;
    
    trial

    init=struct;
    init.ve=sqrt(75)*randn(1,300)+dm(trial-1);
    init.vi=sqrt(80)*randn(1,100);
    m=max(max(init.ve),max(init.vi))-100;
    init.ve=init.ve-m-1;
    init.vi=init.vi-m-1;
    init.he=zeros(2,300);
    init.hi=zeros(2,100);

    res_lif_small=model_LIF_FC(param2,init);
    ind=find(res_lif_small.wave_spike_count>300);
    ind=min(ind);
    if res_lif_small.MFE_time(ind,1)>10
        pos=ceil(res_lif_small.MFE_time(ind,1)*10-20);
        dm(trial)=mean(res_lif_small.VE(pos,:))-mean(res_lif_small.VI(pos,:));
    else
        pos=ceil(res_lif_small.MFE_time(ind+1,1)*10-20);
        dm(trial)=mean(res_lif_small.VE(pos,:))-mean(res_lif_small.VI(pos,:));
    end
    res{trial}=res_lif_small;
    ind(trial)=pos;
end


subplot(1,3,1);

scatter(dm(1:length(dm)-1),dm(2:length(dm)),'.');
grid on;
hold on;
plot([-30:30],[-30:30]);
xlim([-15 30]);
ylim([-15 30]);
xticks([-14.99,-10:5:30])
xticklabels({'-0.15','-0.1','-0.05','0','0.05','0.1','0.15','0.2','0.25','0.3'});
yticks([-14.99,-10:5:30])
yticklabels({'-0.15','-0.1','-0.05','0','0.05','0.1','0.15','0.2','0.25','0.3'});
set(gca,'FontSize',8.5)

ylabel('$$\hat{m}_1$$','FontSize',11,'Interpreter','latex');
xlabel('$m_0$','FontSize',11,'Interpreter','latex');

dm=zeros(1,300);
res=cell(1,300);
index=zeros(1,300);



for trial=2:length(dm)
    param2=param;
    param2.duration=0.08;
    param2.s_ei     = 4.91*0.44;
    
    trial

    init=struct;
    init.ve=sqrt(75)*randn(1,300)+dm(trial-1);
    init.vi=sqrt(80)*randn(1,100);
    m=max(max(init.ve),max(init.vi))-100;
    init.ve=init.ve-m-1;
    init.vi=init.vi-m-1;
    init.he=zeros(2,300);
    init.hi=zeros(2,100);

    res_lif_small=model_LIF_FC(param2,init);
    ind=find(res_lif_small.wave_spike_count>300);
    ind=min(ind);
    if res_lif_small.MFE_time(ind,1)>10
        pos=ceil(res_lif_small.MFE_time(ind,1)*10-20);
        dm(trial)=mean(res_lif_small.VE(pos,:))-mean(res_lif_small.VI(pos,:));
    else
        pos=ceil(res_lif_small.MFE_time(ind+1,1)*10-20);
        dm(trial)=mean(res_lif_small.VE(pos,:))-mean(res_lif_small.VI(pos,:));
    end
    res{trial}=res_lif_small;
    ind(trial)=pos;
end


subplot(1,3,2);

scatter(dm(1:length(dm)-1),dm(2:length(dm)),'.');
grid on;
hold on;
plot([-30:30],[-30:30]);
xlim([-15 30]);
ylim([-15 30]);
xticks([-14.99,-10:5:30])
xticklabels({'-0.15','-0.1','-0.05','0','0.05','0.1','0.15','0.2','0.25','0.3'});
yticks([-14.99,-10:5:30])
yticklabels({'-0.15','-0.1','-0.05','0','0.05','0.1','0.15','0.2','0.25','0.3'});
set(gca,'FontSize',8.5)

ylabel('$$\hat{m}_1$$','FontSize',11,'Interpreter','latex');
xlabel('$m_0$','FontSize',11,'Interpreter','latex');

dm=zeros(1,300);
res=cell(1,300);
index=zeros(1,300);



for trial=2:length(dm)
    param2=param;
    param2.duration=0.08;
    param2.s_ei     = 4.91*0.421;
    
    trial

    init=struct;
    init.ve=sqrt(75)*randn(1,300)+dm(trial-1);
    init.vi=sqrt(80)*randn(1,100);
    m=max(max(init.ve),max(init.vi))-100;
    init.ve=init.ve-m-1;
    init.vi=init.vi-m-1;
    init.he=zeros(2,300);
    init.hi=zeros(2,100);

    res_lif_small=model_LIF_FC(param2,init);
    ind=find(res_lif_small.wave_spike_count>300);
    ind=min(ind);
    if res_lif_small.MFE_time(ind,1)>10
        pos=ceil(res_lif_small.MFE_time(ind,1)*10-20);
        dm(trial)=mean(res_lif_small.VE(pos,:))-mean(res_lif_small.VI(pos,:));
    else
        pos=ceil(res_lif_small.MFE_time(ind+1,1)*10-20);
        dm(trial)=mean(res_lif_small.VE(pos,:))-mean(res_lif_small.VI(pos,:));
    end
    res{trial}=res_lif_small;
    ind(trial)=pos;
end


subplot(1,3,3);

scatter(dm(1:length(dm)-1),dm(2:length(dm)),'.');
grid on;
hold on;
plot([-30:30],[-30:30]);
xlim([-15 30]);
ylim([-15 30]);
xticks([-14.99,-10:5:30])
xticklabels({'-0.15','-0.1','-0.05','0','0.05','0.1','0.15','0.2','0.25','0.3'});
yticks([-14.99,-10:5:30])
yticklabels({'-0.15','-0.1','-0.05','0','0.05','0.1','0.15','0.2','0.25','0.3'});
set(gca,'FontSize',8.5)

ylabel('$$\hat{m}_1$$','FontSize',11,'Interpreter','latex');
xlabel('$m_0$','FontSize',11,'Interpreter','latex');



set(gcf,'Position',[10,10,1200,400]);

% saveas(gcf ,['outputs/dm mapping graph/parameter_searching/',name,'.fig']);
% close(gcf);
end

%%