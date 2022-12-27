%% Setting Path

addpath('module');

%% Parameters for Model LIF

alpha = 1;
beta = 1;
param.ne = 300;
param.ni = 100;
param.M        = 100*beta;
param.Mr       = 66*beta;
param.p_ee     = 1;
param.p_ie     = 1;
param.p_ei     = 1;
param.p_ii     = 1;
param.s_ee     = 5*0.15/alpha*beta;
param.s_ie     = 2*0.5/alpha*beta;
param.s_ei     = 4.91*0.424/alpha*beta;
param.s_ii     = 4.91*0.40/alpha*beta;

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

% Model LIF

tic;
res_lif=model_LIF3(param,[]);
toc;
figure;
rasterplot(res_lif, param);
xlim([2000, 3000]);

%% MFE figure

figure;
rasterplot(res_lif, param);
hold on;
xline(res_lif.MFE_time(1:res_lif.wave_count,1),'b','LineWidth',1);
hold on;
xline(res_lif.MFE_time(1:res_lif.wave_count,2),'LineWidth',1);
xlim([2000, 3000]);

%% Period 1
param.duration = 30;
param.s_ei = 4.91*0.405/alpha*beta;
tic;
res_lif1=model_LIF3(param,[]);
toc;


%% Period 2

param.s_ei = 4.91*0.44/alpha*beta;
tic;
res_lif2=model_LIF3(param,[]);
toc;


%% Period 3
beta = 3;

param.M        = 100*beta;
param.Mr       = 66*beta;
param.s_ee     = 5*0.15/alpha*beta;
param.s_ie     = 2*0.5/alpha*beta;
param.s_ei     = 4.91*0.418/alpha*beta;
param.s_ii     = 4.91*0.40/alpha*beta;

param.lambda_e = 7000*beta;
param.lambda_i = 7000*beta;
tic;
res_lif3=model_LIF3(param,[]);
toc;

%% Raster
figure;
subplot(3,1,1);
rasterplot(res_lif1, param);
xlabel('Time (ms)');
xlim([2500,3000]);
set(gca,'fontsize',15,'fontname','Arial');
subplot(3,1,2);
rasterplot(res_lif2, param);
xlabel('Time (ms)');
xlim([2500,3000]);
set(gca,'fontsize',15,'fontname','Arial');
subplot(3,1,3);
rasterplot(res_lif3, param);
xlabel('Time (ms)');
xlim([2500,3000]);
set(gcf,'Position',[10,10,1200,900]);
set(gca,'fontsize',15,'fontname','Arial');
%print(gcf, '-dpdf', 'figure/Publication/figure1-b-1.pdf','-bestfit');

%% Save Data
save('fig1','res_lif1', 'res_lif2','res_lif3','sd1','sd2','sd3','spec1','spec2','spec3')

%% Spec
param.sdbin = 2.5;
param.spectrogram_timewindow = 200;
param.frequency_range = [10,80];

sd1 = spikedensity(res_lif1, param);
sd2 = spikedensity(res_lif2, param);
sd3 = spikedensity(res_lif3, param);

figure;
subplot(3,1,1);
spec1 = spectrogram(sd1.e, param);
subplot(3,1,2);
spec2 = spectrogram(sd2.e, param);
subplot(3,1,3);
spec3 = spectrogram(sd3.e, param);
set(gcf,'Position',[10,10,1200,900]);


%% Spec Plot
m_spec1 = mean(spec1, 2);
m_spec2 = mean(spec2, 2);
m_spec3 = mean(spec3, 2);
d_spec1 = zeros(71, 10);
d_spec2 = zeros(71, 10);
d_spec3 = zeros(71, 10);
for i = 1:9
d_spec1(:,i) = mean(spec1(:,(i-1)*1192+1:i*1192),2);
d_spec2(:,i) = mean(spec2(:,(i-1)*1192+1:i*1192),2);
d_spec3(:,i) = mean(spec3(:,(i-1)*1192+1:i*1192),2);
end 
d_spec1(:,10) = mean(spec1(:,9*1192+1:11922),2);
d_spec2(:,10) = mean(spec2(:,9*1192+1:11922),2);
d_spec3(:,10) = mean(spec3(:,9*1192+1:11922),2);
var_spec1 = sum((d_spec1 - mean(d_spec1,2)).^2,2);
var_spec2 = sum((d_spec2 - mean(d_spec2,2)).^2,2);
var_spec3 = sum((d_spec3 - mean(d_spec3,2)).^2,2);
se_spec1 = sqrt(var_spec1/10);
se_spec2 = sqrt(var_spec2/10);
se_spec3 = sqrt(var_spec3/10);

%%
fre  = 10:1:80;
figure;
errorbar(fre, m_spec1, se_spec1);
hold on;
set(gca,'fontsize',15,'fontname','Arial');
errorbar(fre, m_spec2, se_spec2);
hold on;
set(gca,'fontsize',15,'fontname','Arial');
errorbar(fre, m_spec3, se_spec3);
xlim([10,70]);
xlabel('Freq (Hz)')
ylabel('Avg Spec')
legend('Period 1', 'Period 2', 'Period 3');
set(gcf,'Position',[10,10,1200,300]);
set(gca,'fontsize',15,'fontname','Arial');
%print(gcf, '-dpdf', 'figure/Publication/figure1-b-3.pdf','-bestfit');
%% Figure 2

%% a
beta = 1;

param.M        = 100*beta;
param.Mr       = 66*beta;
param.lambda_e = 7000*beta;
param.lambda_i = 7000*beta;
param.p_ee     = 0.15;
param.p_ie     = 0.5;
param.p_ei     = 0.42;
param.p_ii     = 0.4;
param.s_ee     = 5;
param.s_ie     = 2;
param.s_ei     = 4.91;
param.s_ii     = 4.91;
param.tau_ee   = 1.4;
param.tau_ie   = 1.2;
param.duration = 3;
tic;
res_lif=model_LIF3(param,[]);
toc;

beta = 1;

param.M        = 100*beta;
param.Mr       = 66*beta;
param.lambda_e = 7000*beta;
param.lambda_i = 7000*beta;
param.p_ee     = 1;
param.p_ie     = 1;
param.p_ei     = 1;
param.p_ii     = 1;
param.s_ee     = 5*0.15;
param.s_ie     = 2*0.5;
param.s_ii     = 4.91*0.40;
param.s_ei = 4.91*0.43;
tic;
res_lif4=model_LIF3(param,[]);
toc;
%%
figure;
subplot(2,1,1);
rasterplot(res_lif, param);
hold on;
xline(res_lif.MFE_time(1:res_lif.wave_count,1),'b','LineWidth',1);
hold on;
xline(res_lif.MFE_time(1:res_lif.wave_count,2),'LineWidth',1);
xlim([2500, 3000]);
xlabel('Time (ms)')
set(gca,'fontsize',12,'fontname','Arial');

subplot(2,1,2);
rasterplot(res_lif4, param);
hold on;
xline(res_lif4.MFE_time(1:res_lif4.wave_count,1),'b','LineWidth',1);
hold on;
xline(res_lif4.MFE_time(1:res_lif4.wave_count,2),'LineWidth',1);
xlim([2500, 3000]);
xlabel('Time (ms)')
set(gca,'fontsize',12,'fontname','Arial');
set(gcf,'Position',[10,10,1200,600]);
%print(gcf, '-dpdf', 'figure/Publication/figure2-1.pdf','-bestfit');

