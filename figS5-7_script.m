%% Setting Path
addpath('module');
%% Parameters for Model LIF

alpha = 1;
beta = 3;
param.ne = 300;
param.ni = 100;
param.M        = 100;
param.Mr       = 66;
% param.p_ee     = 0.15;
% param.p_ie     = 0.5;
% param.p_ei     = 0.415;
% param.p_ii     = 0.4;
param.p_ee     = 0.3;
param.p_ie     = 0.3;
param.p_ei     = 0.3;
param.p_ii     = 0.3;
% param.s_ee     = 5/alpha*beta;
% param.s_ie     = 2/alpha*beta;
% param.s_ei     = 4.91/alpha*beta;
% param.s_ii     = 4.91/alpha*beta;
param.s_ee     = 5*0.15/alpha/param.p_ee;
param.s_ie     = 2*0.5/alpha/param.p_ie;
param.s_ei     = 4.91*0.416  /alpha/param.p_ei;
param.s_ii     = 4.91*0.4/alpha/param.p_ii;
param.s_exe = 1/beta;
param.s_exi = 1/beta;

param.lambda_e = 7000*beta;
param.lambda_i = 7000*beta;
param.tau_ee   = 1.4;
param.tau_ie   = 1.2;
param.tau_ei    = 4.5;
param.tau_ii    = 4.5;
param.tau_re=0;
param.tau_ri=0;
param.duration = 10;

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
figure;
tic;
res_lif1=model_LIF3(param,[]);
toc;
figure;
rasterplot(res_lif1, param);
xlim([2000, 3000]);


%%
tic;
for j=1:5
    res_lif1=model_LIF3_FC(param,[]);
    mt=res_lif1.MFE_time;
    s=find(mt(:,1)>1000,1);
    p=find(mt(:,1)==0,1);

    subplot(5,3,(j-1)*3+1);
    histogram(mt(s+1:p-1,1)-mt(s:p-2,2),[0:0.25:25],'Normalization','probability');
    if j==1
        title('Length of Inter-MFE Intervals')
        
    end

    if j==5
        xlabel('Time (ms)')
    end
    ylim([0,0.2]);
    if j==3
    ylabel('Empirical Probability',FontSize=17);
    end
    set(gca,'FontSize',15);

    subplot(5,3,(j-1)*3+2);
    histogram(mt(s:p-2,2)-mt(s:p-2,1),[0:0.25:15],'Normalization','probability');
    if j==1
        title('Length of MFE');
    end
    if j==5
        xlabel('Time (ms)')
    end

    if j==3
    ylabel('Empirical Probability',FontSize=17);
    end
    ylim([0 0.2])
    set(gca,'FontSize',15);
    delta_time = 5;
    delta_time = delta_time / 1000;
    size1 = floor(param.duration/delta_time);
    t = 1:1:size1;
    t = t .* delta_time .* 1000 - delta_time/2;
    sd_lif = zeros(size1,1);
    spike_E = res_lif1.spike(:,1:300);
    spike_I = res_lif1.spike(:,301:400);
    FrE=fr_rate(spike_E,[1,0.5,10]);
    FrI=fr_rate(spike_I,[1,0.5,10]);
    spike_E = spike_E(2:max(spike_E(1,:))+1,:);
    for i = 1:size1
        sd_lif(i) = sum(sum(spike_E<= i*delta_time & spike_E> (i-1)*delta_time));
    end


    delta_time=5;
    amp=zeros(p-1-s,1);
    for i=s:p-2
        amp(i-s+1)=sum(sd_lif(ceil(mt(i,1)/delta_time):ceil(mt(i,2)/delta_time)));
    end

    
    subplot(5,3,(j-1)*3+3);
    histogram(amp,[0:5:500],'Normalization','probability');
    if j==1
        title('Amplitute of MFE');
    end
    if j==5
        xlabel('Spikecount')
    end
    if j==3
    ylabel('Empirical Probability',FontSize=17);
    end
    ylim([0 0.3]);
    yticks([0 0.1 0.2 0.3]);
    set(gca,'FontSize',15);

%     subplot(5,4,(j-1)*4+1);
%     hold on
%     color_matrix = [1,0.38,0.27;0,0.75,1];
%     b = bar(1,mean(FrE),0.75,'stacked');
%     set(b(1),'facecolor',color_matrix(1,:))
%     er=errorbar(1,mean(FrE),sqrt(var(FrE))/sqrt(length(FrE)));
%     er.Color = [0 0 0];                            
%     er.LineStyle = 'none';
%     er.CapSize = 18;
%     er.LineWidth = 1;
% 
%     b = bar(2,mean(FrI),0.75,'stacked');
%     set(b(1),'facecolor',color_matrix(2,:))
%     er=errorbar(2,mean(FrI),sqrt(var(FrI))/sqrt(length(FrI)));
%     er.Color = [0 0 0];                            
%     er.LineStyle = 'none';
%     er.CapSize = 18;
%     er.LineWidth = 1;
% 
%     box on;
%     set(gca,'XTick',[1 2]);
%     xticklabels({'FrE','FrI'})
%     if j==1
%         title('Firing Rate')
%         ylabel('spikes/sec')
%     end
% 
%     if j==5
%         xlabel('Firing Rate of E/I Neurons (spikes/sec)')
%     end

disp(['FrE=',num2str(mean(FrE)),'\pm',num2str(sqrt(var(FrE))/sqrt(length(FrE)))]);
disp(['FrI=',num2str(mean(FrI)),'\pm',num2str(sqrt(var(FrI))/sqrt(length(FrI)))]);




end
set(gcf,'Position',[10,10,1500,2000]);

toc;


%% 原LIFmodel图
tic;
for j=1
    res_lif1=model_LIF3(param,[]);
    mt=res_lif1.MFE_time;
    s=find(mt(:,1)>1000,1);
    p=find(mt(:,1)==0,1);

    subplot(1,3,(j-1)*3+1);
    histogram(mt(s+1:p-1,1)-mt(s:p-2,2),[0:0.25:25],'Normalization','probability');
    if j==1
        title('Length of Inter-MFE Intervals')
        
    end

    if j==1
        xlabel('Time (ms)')
    end
    ylim([0,0.2]);
    if j==1
    ylabel('Empirical Probability',FontSize=17);
    end
    set(gca,'FontSize',15);

    subplot(1,3,(j-1)*3+2);
    histogram(mt(s:p-2,2)-mt(s:p-2,1),[0:0.25:15],'Normalization','probability');
    if j==1
        title('Length of MFE');
    end
    if j==1
        xlabel('Time (ms)')
    end

    if j==1
    ylabel('Empirical Probability',FontSize=17);
    end
    ylim([0 0.2])
    set(gca,'FontSize',15);
    delta_time = 5;
    delta_time = delta_time / 1000;
    size1 = floor(param.duration/delta_time);
    t = 1:1:size1;
    t = t .* delta_time .* 1000 - delta_time/2;
    sd_lif = zeros(size1,1);
    spike_E = res_lif1.spike(:,1:300);
    spike_I = res_lif1.spike(:,301:400);
    FrE=fr_rate(spike_E,[1,0.5,10]);
    FrI=fr_rate(spike_I,[1,0.5,10]);
    spike_E = spike_E(2:max(spike_E(1,:))+1,:);
    for i = 1:size1
        sd_lif(i) = sum(sum(spike_E<= i*delta_time & spike_E> (i-1)*delta_time));
    end


    delta_time=5;
    amp=zeros(p-1-s,1);
    for i=s:p-2
        amp(i-s+1)=sum(sd_lif(ceil(mt(i,1)/delta_time):ceil(mt(i,2)/delta_time)));
    end

    
    subplot(1,3,(j-1)*3+3);
    histogram(amp,[0:5:500],'Normalization','probability');
    if j==1
        title('Amplitute of MFE');
    end
    if j==1
        xlabel('Spikecount')
    end
    if j==1
    ylabel('Empirical Probability',FontSize=17);
    end
    ylim([0 0.3]);
    yticks([0 0.1 0.2 0.3]);
    set(gca,'FontSize',15);

%     subplot(5,4,(j-1)*4+1);
%     hold on
%     color_matrix = [1,0.38,0.27;0,0.75,1];
%     b = bar(1,mean(FrE),0.75,'stacked');
%     set(b(1),'facecolor',color_matrix(1,:))
%     er=errorbar(1,mean(FrE),sqrt(var(FrE))/sqrt(length(FrE)));
%     er.Color = [0 0 0];                            
%     er.LineStyle = 'none';
%     er.CapSize = 18;
%     er.LineWidth = 1;
% 
%     b = bar(2,mean(FrI),0.75,'stacked');
%     set(b(1),'facecolor',color_matrix(2,:))
%     er=errorbar(2,mean(FrI),sqrt(var(FrI))/sqrt(length(FrI)));
%     er.Color = [0 0 0];                            
%     er.LineStyle = 'none';
%     er.CapSize = 18;
%     er.LineWidth = 1;
% 
%     box on;
%     set(gca,'XTick',[1 2]);
%     xticklabels({'FrE','FrI'})
%     if j==1
%         title('Firing Rate')
%         ylabel('spikes/sec')
%     end
% 
%     if j==5
%         xlabel('Firing Rate of E/I Neurons (spikes/sec)')
%     end

disp(['FrE=',num2str(mean(FrE)),'\pm',num2str(sqrt(var(FrE))/sqrt(length(FrE)))]);
disp(['FrI=',num2str(mean(FrI)),'\pm',num2str(sqrt(var(FrI))/sqrt(length(FrI)))]);




end
set(gcf,'Position',[10,10,1500,250]);

toc;

%%
FrE=fr_rate(spike_E,[1,0.5,10]);
sqrt(var(FrE))/sqrt(length(FrE))

%%
function fr=fr_rate(spike,param) %param=[start,end,bin]
    
    start=param(1);
    nd=param(3);
    bin=param(2);
    size=round((nd-start)/bin);
    fr=zeros(1,size);
    for i=1:size
        ind=find(spike>start+(i-1)*bin & spike<start+i*bin);
        fr(i)=length(ind)/bin/300;
    end

end





























