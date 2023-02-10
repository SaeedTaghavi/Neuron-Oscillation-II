% Fig 1
% Matlab script to generate fig 1. 
% Fig 1 has 3 pannels.
% - Panel 1: the topological graph for models 
% - Panel 2: Rasterplots & Spikecounts to show 3 regimes 
% - Panel 3: The average spectrogram of 3 regimes 

%% Load Modules & General Parameters

addpath('module');
load("param");

%% Generate Params 

param.duration = 30;
p = 0.8;
param.p_ee = p;
param.p_ei = p;
param.p_ii = p;
param.p_ie = p;
alpha = 1;


% Period 1
beta = 3;
param1 = param;
param1.M        = 100*beta;
param1.Mr       = 66*beta;
param1.s_ee     = 5*0.15/alpha*beta/p;
param1.s_ie     = 2*0.5/alpha*beta/p;
param1.s_ei     = 4.91*0.40/alpha*beta/p;
param1.s_ii     = 4.91*0.40/alpha*beta/p;
param1.lambda_e = 7000*beta;
param1.lambda_i = 7000*beta;

% Period 2 
beta = 3;
param2 = param;
param2.M        = 100*beta;
param2.Mr       = 66*beta;
param2.s_ee     = 5*0.15/alpha*beta/p;
param2.s_ie     = 2*0.5/alpha*beta/p;
param2.s_ei     = 4.91*0.427/alpha*beta/p;
param2.s_ii     = 4.91*0.40/alpha*beta/p;
param2.lambda_e = 7000*beta;
param2.lambda_i = 7000*beta;

% Period 3 
beta = 3;
param3 = param;
param3.M        = 100*beta;
param3.Mr       = 66*beta;
param3.s_ee     = 5*0.15/alpha*beta/p;
param3.s_ie     = 2*0.5/alpha*beta/p;
param3.s_ei     = 4.91*0.416/alpha*beta/p;
param3.s_ii     = 4.91*0.40/alpha*beta/p;
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

%% Rasterplots

figure;
subplot(3,1,1);
rasterplot(res_lif1, param1);
xlabel('Time (ms)');
xlim([2000,3000]);
set(gca,'fontsize',15,'fontname','Arial');
subplot(3,1,2);
rasterplot(res_lif2, param2);
xlabel('Time (ms)');
xlim([2000,3000]);
set(gca,'fontsize',15,'fontname','Arial');
subplot(3,1,3);
rasterplot(res_lif3, param3);
xlabel('Time (ms)');
xlim([2000,3000]);
set(gcf,'Position',[10,10,1200,900]);
set(gca,'fontsize',15,'fontname','Arial');
% print(gcf, '-dpdf', 'figure/Publication/figure1-b-1.pdf','-bestfit');

%% Spikecount

delta_time = 5;
delta_time = delta_time / 1000;
size1 = floor(param.duration/delta_time);
t = 1:1:size1;
t = t .* delta_time .* 1000 - delta_time/2;
sd_lif1 = zeros(size1,2);
spike_E = res_lif1.spike(:,1:300);
spike_I = res_lif1.spike(:,301:400);
spike_E = spike_E(1:max(spike_E(1,:))+1,:);
spike_I = spike_I(1:max(spike_I(1,:))+1,:);
for i = 1:size1
    sd_lif1(i,1) = sum(sum(spike_E< i*delta_time & spike_E>= (i-1)*delta_time));
    sd_lif1(i,2) = sum(sum(spike_I< i*delta_time & spike_I>= (i-1)*delta_time));
end

sd_lif2 = zeros(size1,2);
spike_E = res_lif2.spike(:,1:300);
spike_I = res_lif2.spike(:,301:400);
spike_E = spike_E(1:max(spike_E(1,:))+1,:);
spike_I = spike_I(1:max(spike_I(1,:))+1,:);
for i = 1:size1
    sd_lif2(i,1) = sum(sum(spike_E< i*delta_time & spike_E>= (i-1)*delta_time));
    sd_lif2(i,2) = sum(sum(spike_I< i*delta_time & spike_I>= (i-1)*delta_time));
end

sd_lif3 = zeros(size1,2);
spike_E = res_lif3.spike(:,1:300);
spike_I = res_lif3.spike(:,301:400);
spike_E = spike_E(1:max(spike_E(1,:))+1,:);
spike_I = spike_I(1:max(spike_I(1,:))+1,:);
for i = 1:size1
    sd_lif3(i,1) = sum(sum(spike_E< i*delta_time & spike_E>= (i-1)*delta_time));
    sd_lif3(i,2) = sum(sum(spike_I< i*delta_time & spike_I>= (i-1)*delta_time));
end

%% Panel B: Rasterplots & Spikecount

figure;
subplot(3,2,1);
rasterplot(res_lif1, param1);
box on;
xlabel('Time (ms)');
xlim([2500,3000]);
set(gca,'fontsize',15,'fontname','Arial');

subplot(3,2,2);
plot(t, sd_lif1(:,1));
ylabel('Spikecount');
xlabel('Time (ms)');
xlim([2500,3000]);
ylim([0,400]);
set(gca,'fontsize',15,'fontname','Arial');

subplot(3,2,3);
rasterplot(res_lif2, param2);
box on;
xlabel('Time (ms)');
xlim([2500,3000]);
set(gca,'fontsize',15,'fontname','Arial');

subplot(3,2,4);
plot(t, sd_lif2(:,1));
xlabel('Time (ms)');
ylabel('Spikecount');
xlim([2500,3000]);
ylim([0,400]);
set(gca,'fontsize',15,'fontname','Arial');

subplot(3,2,5);
rasterplot(res_lif3, param3);
box on;
xlabel('Time (ms)');
xlim([2500,3000]);
set(gca,'fontsize',15,'fontname','Arial');

subplot(3,2,6);
plot(t, sd_lif3(:,1));
xlabel('Time (ms)');
ylabel('Spikecount');
xlim([2500,3000]);
ylim([0,400]);
set(gca,'fontsize',15,'fontname','Arial');
% sgtitle(['p=',num2str(p),',beta=', num2str(beta)]);
set(gcf,'Position',[10,10,1200,900]);
set(gca,'fontsize',15,'fontname','Arial');
print(gcf, '-dpdf', 'figs/Fig1/figure1-b.pdf','-bestfit');


%% Spectrogram of Firing events

param.sdbin = 2.5;
param.spectrogram_timewindow = 200;
param.frequency_range = [5,80];
param.grid = 1;

sd1 = spikedensity(res_lif1, param);
sd2 = spikedensity(res_lif2, param);
sd3 = spikedensity(res_lif3, param);


spec1 = spectrogram(sd1.e, param);
spec2 = spectrogram(sd2.e, param);
spec3 = spectrogram(sd3.e, param);

num_sample = size(spec1, 2);

%
m_spec1 = mean(spec1, 2);
m_spec2 = mean(spec2, 2);
m_spec3 = mean(spec3, 2);

avg_mat = kron(eye(1192),ones(10,1)/10);
q = zeros(1,1192);
avg_mat = [avg_mat; q; q];

s_spec1 = spec1 * avg_mat;
s_spec2 = spec2 * avg_mat;
s_spec3 = spec3 * avg_mat; 
num_s = size(s_spec1,2);


var_spec1 = sum((s_spec1 - m_spec1).^2,2);
var_spec2 = sum((s_spec2 - m_spec2).^2,2);
var_spec3 = sum((s_spec3 - m_spec3).^2,2);
se_spec1 = sqrt(var_spec1/num_s^2);
se_spec2 = sqrt(var_spec2/num_s^2);
se_spec3 = sqrt(var_spec3/num_s^2);

sm_spec1 = conv(m_spec1, ones(5,1)/5, "same");
sm_spec2 = conv(m_spec2, ones(5,1)/5, "same");
sm_spec3 = conv(m_spec3, ones(5,1)/5, "same");

%% Spec of single neuron potential 


param.sdbin = 1;
param.spectrogram_timewindow = 200;
param.frequency_range = [10,80];
param.grid = 1;
num_ve = 300000;
gg = 10;
ve1 = res_lif1.VE(:,2);
ve1 = ve1(2:end);
ve1 = reshape(ve1, gg, num_ve/gg);
ve1 = mean(ve1,1);
ve1 = [0, ve1];
ve2 = res_lif2.VE(:,1);
ve2 = ve2(2:end);
ve2 = reshape(ve2, gg, num_ve/gg);
ve2 = mean(ve2,1);
ve2 = [0, ve2];
ve3 = res_lif3.VE(:,1);
ve3 = ve3(2:end);
ve3 = reshape(ve3, gg, num_ve/gg);
ve3 = mean(ve3,1);
ve1 = ve1 - mean(ve1);
ve2 = ve2 - mean(ve2);
ve3 = ve3- mean(ve3);
ve3 = [0, ve3];

spec1 = spectrogram(ve1, param);
spec2 = spectrogram(ve2, param);
spec3 = spectrogram(ve3, param);

num_sample = size(spec1, 2);

m_spec1 = mean(spec1, 2);
m_spec2 = mean(spec2, 2);
m_spec3 = mean(spec3, 2);
% m_spec1 = conv(m_spec1, ones(3,1)/3, "same");
% m_spec2 = conv(m_spec2, ones(3,1)/3, "same");
% m_spec3 = conv(m_spec3, ones(3,1)/3, "same");

var_spec1 = sum((spec1 - m_spec1).^2,2);
var_spec2 = sum((spec2 - m_spec2).^2,2);
var_spec3 = sum((spec3 - m_spec3).^2,2);
se_spec1 = sqrt(var_spec1/num_sample^2);
se_spec2 = sqrt(var_spec2/num_sample^2);
se_spec3 = sqrt(var_spec3/num_sample^2);

%% Panel C: Spectrogram
fre  = param.frequency_range(1):1:param.frequency_range(2);

figure;
errorbar(fre, sm_spec1, se_spec1);
hold on;
errorbar(fre, sm_spec2, se_spec2);
hold on;
errorbar(fre, sm_spec3, se_spec3);

xlim([11,57]);
ylim([2,400]);
xlabel('Freq (Hz)');
ylabel('Avg Spec ((spikes/sec)^2/Hz)');
legend('Period 1', 'Period 2', 'Period 3');
set(gcf,'Position',[10,10,1200,300]);
% title(['p=',num2str(p),',beta=', num2str(beta)]);
set(gca,'fontsize',15,'fontname','Arial','YScale','log');
% set(gca,'fontsize',15,'fontname','Arial');
%print(gcf, '-dpdf', 'figure/Publication/figure1-b-3.pdf','-bestfit');

%% Save Data
save('outputs/fig1/data.mat','res_lif1','res_lif2','res_lif3','-v7.3')
