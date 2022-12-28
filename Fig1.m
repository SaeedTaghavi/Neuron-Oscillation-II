% Fig 1
% Matlab script to generate fig 1. 
% Fig 1 has 3 pannels.
% - Panel 1: the topological graph for models 
% - Panel 2: Rasterplots to show 3 regimes 
% - Panel 3: The average spectrogram of 3 regimes 

%% Load Modules & General Parameters

addpath('module');
load("param");

%% Generate Params 

param.duration = 30;
param.p_ee = 0.8;
param.p_ei = 0.8;
param.p_ii = 0.8;
param.p_ie = 0.8;

% Period 1
param1 = param;
param1.s_ei = 4.91*0.405/alpha*beta/0.8;
param1.s_ii = param.s_ii/0.8;
param1.s_ee = param.s_ee/0.8;
param1.s_ie = param.s_ie/0.8;

% Period 2 
param2 = param;
param2.s_ei = 4.91*0.44/alpha*beta/0.8;
param2.s_ii = param.s_ii/0.8;
param2.s_ee = param.s_ee/0.8;
param2.s_ie = param.s_ie/0.8;

% Period 3 
beta = 3;
param3 = param;
param3.M        = 100*beta;
param3.Mr       = 66*beta;
param3.s_ee     = 5*0.15/alpha*beta/0.8;
param3.s_ie     = 2*0.5/alpha*beta/0.8;
param3.s_ei     = 4.91*0.418/alpha*beta/0.8;
param3.s_ii     = 4.91*0.40/alpha*beta/0.8;
param3.lambda_e = 7000*beta;
param3.lambda_i = 7000*beta;

%% Models 

tic;
res_lif1=model_LIF_FC(param1,[]);
toc;

tic;
res_lif2=model_LIF_FC(param2,[]);
toc;

tic;
res_lif3=model_LIF_FC(param3,[]);
toc;

%% Raster
figure;
subplot(3,1,1);
rasterplot(res_lif1, param1);
xlabel('Time (ms)');
xlim([2500,3000]);
set(gca,'fontsize',15,'fontname','Arial');
subplot(3,1,2);
rasterplot(res_lif2, param2);
xlabel('Time (ms)');
xlim([2500,3000]);
set(gca,'fontsize',15,'fontname','Arial');
subplot(3,1,3);
rasterplot(res_lif3, param3);
xlabel('Time (ms)');
xlim([2500,3000]);
set(gcf,'Position',[10,10,1200,900]);
set(gca,'fontsize',15,'fontname','Arial');
% print(gcf, '-dpdf', 'figure/Publication/figure1-b-1.pdf','-bestfit');

%% Save Data
save('output/fig1','res_lif1', 'res_lif2','res_lif3','sd1','sd2','sd3','spec1','spec2','spec3')

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
