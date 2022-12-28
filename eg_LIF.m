%% Setting Path

addpath('module');

%% Parameters for Model LIF

param.ne = 300;
param.ni = 100;
param.M        = 100;
param.Mr       = 66;
param.p_ee     = 0.15;
param.p_ie     = 0.5;
param.p_ei     = 0.5;
param.p_ii     = 0.4;
% param.p_ee     = 1;
% param.p_ie     = 1;
% param.p_ei     = 1;
% param.p_ii     = 1;
param.s_ee     = 5;
param.s_ie     = 2;
param.s_ei     = 4.91;
param.s_ii     = 4.91;
% param.s_ee     = 5*0.15;
% param.s_ie     = 2*0.5;
% param.s_ei     = 4.91*0.45;
% param.s_ii     = 4.91*0.4;

param.lambda_e = 7000;
param.lambda_i = 7000;
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

param2=param;
param2.gridsize=0.1;

%% Model LIF

tic;
res_lif=model_LIF3(param2);
toc;

%%

spikecount_e=res_lif.spikecount_e;
spikecount_i=res_lif.spikecount_i;
mfe_init=[];
i=1;
while i<length(spikecount_e)
    if spikecount_e(i)~=0
        mfe_init=[mfe_init,i];
        i=i+50;
    else
    	i=i+1;
    end
end

dm=zeros(1,length(mfe_init));
for i=1:length(mfe_init)
    ind=mfe_init(i)-3;
    dm(i)=mean(res_lif.VE(ind,:))-mean(res_lif.VI(ind,:));
end

%%
figure;
subplot(3,4,[1,2,5,6,9,10]);
scatter(dm(1:end-1),dm(2:end));
xlabel('m_{VE} - m_{VI}');
ylabel('m_{VE} - m_{VI}');
subplot(3,4,[3,4]);
yyaxis left;
plot((mfe_init-3)*param2.gridsize,dm(1:end));
ylabel('m_{VE} - m_{VI}');
yyaxis right;
t = 0:param2.gridsize:param2.duration*1000;
plot(t, res_lif.spikecount_e);
ylabel('Spikecount');
xlabel('Time (ms)');
xlim([2000,3000]);
subplot(3,4,[7,8]);
dmi = dm>10;
plot((mfe_init-3)*param2.gridsize,dmi);
xlim([2000,3000]);
xlabel('Time (ms)');
ylabel('m_{VE}-m_{VI}>10');
sgtitle('LIF P_{EI} = 0.421');
subplot(3,4,[11,12]);
rasterplot(res_lif,param2);
xlabel('Time (ms)');
xlim([2000,3000]);
set(gcf,'Position',[10,10,1200,600]);

%% Single Experiment
PEI = 0.43;
param2.s_ei = PEI*4.91;

tic;
res_lif=model_LIF(param2);
toc;

spikecount_e=res_lif.spikecount_e;
spikecount_i=res_lif.spikecount_i;
mfe_init=[];
i=1;
while i<length(spikecount_i)
    if spikecount_i(i)~=0
        mfe_init=[mfe_init,i];
        i=i+500;
    else
    	i=i+1;
    end
end

dm=zeros(1,length(mfe_init));
for i=1:length(mfe_init)
    ind=mfe_init(i)-3;
    dm(i)=mean(res_lif.VE(ind,:))-mean(res_lif.VI(ind,:));
end

dma = mean(res_lif.VE,2) - mean(res_lif.VI,2);
t = 0:param2.gridsize:param2.duration*1000;

%%

figure;
subplot(2,4,[1,2,5,6]);
scatter(dm(1:end-1),dm(2:end));
xlabel('m_{VE} - m_{VI}');
ylabel('m_{VE} - m_{VI}');
subplot(2,4,[3,4]);
yyaxis left;
plot(t,dma);
ylabel('m_{VE} - m_{VI}');
yyaxis right;
t = 0:param2.gridsize:param2.duration*1000;
plot(t, res_lif.spikecount_e);
ylabel('Spikecount');
xlabel('Time (ms)');
xlim([2000,3000]);
subplot(2,4,[7,8]);
rasterplot(res_lif,param2);
xlabel('Time (ms)');
xlim([2000,3000]);
set(gcf,'Position',[10,10,1200,600]);
name1 = ['LIF-P_{EI} = ', num2str(PEI)];
name2 = ['LIF-p_ei=',num2str(PEI)];
sgtitle(name1);

saveas(gcf,['figure/LIF/',name2,'.jpg']);
    saveas(gcf,['figure/LIF/',name2,'.fig']);