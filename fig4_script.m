%% Setting Path

addpath('module');

%% Parameters

ps = 0.25:0.03:0.46;
seis = 4.91*(0.40:0.004:0.444);
taus = 0:0.5:5.5;
sexs = 0.2:0.8:9;

ps=0.25:0.25:1
seis=1.96:0.06:2.14
NGEs = cell(12,1);
NGIs = cell(12,1);
HEs = cell(12,1);
HIs = cell(12,1);
RESs = cell(12,1);
MEs = cell(12,1);
MIs = cell(12,1);

% Parameters for Model LIF
for i=1:4



%    p=ps(i);
      p=0.46;
    alpha = 1;
    beta = 1;
    param.ne = 300;
    param.ni = 100;
    param.M        = 100*beta;
    param.Mr       = 66*beta;
    param.p_ee     = p;
    param.p_ie     = p;
    param.p_ei     = p;
    param.p_ii     = p;
    param.s_ee     = 5*0.15/alpha*beta/p;
    param.s_ie     = 2*0.5/alpha*beta/p;
    param.s_ei     = 4.91*0.421/alpha*beta/p;
%    param.s_ei     = seis(i)/alpha*beta/p;
    param.s_ii     = 4.91*0.40/alpha*beta/p;



    param.tau_ee   = 1.4;
    param.tau_ie   = 1.2;
    param.tau_ei    = 4.5;
    param.tau_ii    = 4.5;
%     param.tau_re=0;
%     param.tau_ri=0;
    param.tau_re=taus((i-1)*3+1);
    param.tau_ri=taus((i-1)*3+1);
    param.duration = 3;
    param.delta_time = 0.1;

    param.factor_Leak=0;
    param.LeakE = 20;
    param.s_exe = 1;
    param.s_exi = 1;
    param.lambda_e = 7000*beta/sexs(2);
    param.lambda_i = 7000*beta/sexs(2);
%     param.s_exe = sexs((i-1)*3+1);
%     param.s_exi = sexs((i-1)*3+1);
%     param.lambda_e = 7000*beta/sexs((i-1)*3+1);
%     param.lambda_i = 7000*beta/sexs((i-1)*3+1);
    param.LeakI = 16.7;
    param.factor_Leak = inf;

    param2=param;
    param2.gridsize=0.1;

    % Model LIF

    tic;
    res_lif=model_LIF3(param2,[]);
    toc;

    VE = res_lif.VE;
    VI = res_lif.VI;
    HE = res_lif.HE;
    HI = res_lif.HI;
    HE = sum(HE,2);
    HI = sum(HI,2);
    VE = VE./beta;
    VI = VI./beta;
    NGE = sum(VE>50,2);
    NGI = sum(VI>50,2);
    ME = mean(VE,2);
    MI = mean(VI,2);

    NGEs{i} = NGE;
    NGIs{i} = NGI;
    HEs{i} = HE;
    HIs{i} = HI;
    RESs{i} = res_lif;
    MEs{i} = ME;
    MIs{i} = MI;
    % figure;
    % a=plot3(HE, HI, NGE);
    % a.Color(4) = .2;
    % box on;
    % grid on;
    % xlabel('HE');
    % ylabel('HI');
    % zlabel('NGE');
end
c=turbo;
%%
figure;
for i = 1:12
    subplot(3,4,i);
    mc_plot(RESs{i},HEs{i}, HIs{i}, NGEs{i},c);
    %a.Color(4) = .2;
    box on;
    grid on;
    xlabel('HE');
    ylabel('HI');
    zlabel('NGE');
    %     title(['S_{ex}=',num2str(sexs(i))]);
    %     title(['S_{ei}=',num2str(seis(i))]);
    %     title(['P=',num2str(ps(i))]);
    title(['Ref=',num2str(taus((i-1)*3+1))]);
    view([-30,40]);
end
set(gcf,'Position',[10,10,1600,300]);
%saveas(gcf,['P=', num2str(ps(1)),'-',num2str(ps(12)),'.fig'])

%%
figure;
for i = 1:4
    subplot(1,4,i);
    dm_record=mc_plot(RESs{i},HEs{i}, HIs{i}, NGEs{i},c);
    %a.Color(4) = .2;
    box on;
    grid on;
    xlabel('g^{E}')
    ylabel('g^{I}')
    if i ==1
        zlabel('# of gate E neurons')
    end
%     xticks([0 12000 24000]*i);
%     xticklabels([0 120 240]*i);
%     yticks([0 6000 12000 18000]*i);
%     yticklabels([0 60 120 180]*i);
    xticks([0 20000 40000]);
    xticklabels([0 200 400]);
    yticks([0 10000 20000 30000]);
    yticklabels([0 100 200 300]);
         title(['S^{ext}=',num2str(sexs((i-1)*3+1)),'×','10^{-2}']);

%      title(['S^{EI}=',num2str(seis(i)),'×','10^{-2}']);
%         title(['P=',num2str(ps(i))]);
    view([-25,45]);
    xlim([0 40000]);
    ylim([0 30000]);
%          xlim([0 24000*i]);
%          ylim([0 18000*i]);
    %      xticks([0 12000 24000]*i);
    %      yticks([0 6000 12000 18000]*i)
end
set(gcf,'Position',[10,10,1600,300]);
%saveas(gcf,['P=', num2str(ps(1)),'-',num2str(ps(12)),'.fig'])

%%
figure;
for i = 1:4
    subplot(1,4,i);
    dm_record=start_point_plot(RESs{i},MEs{i},MIs{i},HIs{i},c);
    %a.Color(4) = .2;
    box on;
    grid on;
    xlabel('m^{E}')
    ylabel('m^{I}')
    if i ==1
        zlabel('g^{I}')
    end
%     xticks([0 12000 24000]*i);
%     xticklabels([0 120 240]*i);
%     yticks([0 6000 12000 18000]*i);
%     yticklabels([0 60 120 180]*i);
%     xticks([0 20000 40000]);
%     xticklabels([0 200 400]);
    xticks([50 55 60 65 70 75 80 85 90]);
    yticks([30 40 50 60 70 80 90])

    zticks([0 2000 4000 6000]);
    zticklabels([0 20 40 60]);
%          title(['S^{ext}=',num2str(sexs((i-1)*3+1)),'×','10^{-2}']);
%         title(['S_{ei}=',num2str(seis(i))]);
%      title(['S^{EI}=',num2str(seis(i)),'×','10^{-2}']);
%          title(['P=',num2str(ps(i))]);
   title(['\tau_R=',num2str(taus((i-1)*3+1))])
    view([-25,55]);
%       xlim([50 90]);
%       ylim([50 90]);
      zlim([0 4000]);
%          xlim([0 24000*i]);
%          ylim([0 18000*i]);
    %      xticks([0 12000 24000]*i);
    %      yticks([0 6000 12000 18000]*i)
end
set(gcf,'Position',[10,10,1600,300]);
%saveas(gcf,['P=', num2str(ps(1)),'-',num2str(ps(12)),'.fig'])
%%
p1= 0.13;
p2 = 1;
alpha = 1;
beta = 1;
param.ne = 300;
param.ni = 100;
param.M        = 100*beta;
param.Mr       = 66*beta;
param.p_ee     = p1;
param.p_ie     = p2;
param.p_ei     = p2;
param.p_ii     = p2;

param.s_ee     = 5*0.15/alpha*beta/p1;
param.s_ie     = 2*0.5/alpha*beta/p2;
param.s_ei     = 4.91*0.419/alpha*beta/p2;
param.s_ii     = 4.91*0.40/alpha*beta/p2;

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

param2=param;
param2.gridsize=0.1;

% Model LIF

tic;
res_lif=model_LIF3(param2,[]);
toc;

VE = res_lif.VE;
VI = res_lif.VI;
HE = res_lif.HE;
HI = res_lif.HI;
HE = sum(HE,2);
HI = sum(HI,2);
VE = VE./10;
VI = VI./10;
NGE = sum(VE>60,2);
NGI = sum(VI>60,2);

%
figure;
a=plot3(HE, HI, NGE);
a.Color(4) = .2;
box on;
grid on;
xlabel('HE');
ylabel('HI');
zlabel('NGE');
title(['P=',num2str(p1)]);
view([-30,40]);
%%
i=2
mc_plot(RESs{i},HEs{i}, HIs{i},NGEs{i},c);
%%
function dm_record=mc_plot(res_lif,HE,HI,NGE,c)
ind=find(res_lif.MFE_time(:,1)~=0,1,"last")-1;
for i=1:ind
    if res_lif.wave_spike_count(i)>150
        t_start=round(res_lif.MFE_time(i,1)*10);
        t_end=round(res_lif.MFE_time(i+1,1)*10);
        dm=mean(res_lif.VE(t_start-20,:))-mean(res_lif.VI(t_start-20,:));
        dm=ceil(dm*2.8+150);
        dm_record(i)=dm;
        plot3(HE(t_start:t_end),HI(t_start:t_end),NGE(t_start:t_end),'Color',c(dm,:));
        plot3(HE(t_start),HI(t_start),NGE(t_start),'.','Color','c','MarkerSize',15)
        hold on;
        
    else
        t_start=round(res_lif.MFE_time(i,1)*10);
        t_end=round(res_lif.MFE_time(i+1,1)*10);
        plot3(HE(t_start:t_end),HI(t_start:t_end),NGE(t_start:t_end),'b');
        plot3(HE(t_start),HI(t_start),NGE(t_start),'.','Color','c','MarkerSize',15)
        hold on;
        
    end
end
 zlim([0 300])
end

%%
function dm_record=start_point_plot(res_lif,HE,HI,NGE,c)
ind=find(res_lif.MFE_time(:,1)~=0,1,"last")-1;
for i=1:ind
    if res_lif.wave_spike_count(i)>150
        t_start=round(res_lif.MFE_time(i,1)*10);
        t_end=round(res_lif.MFE_time(i+1,1)*10);
        dm=mean(res_lif.VE(t_start-20,:))-mean(res_lif.VI(t_start-20,:));
        dm=ceil(dm*2.8+150);
        dm_record(i)=dm;
        plot3(HE(t_start),HI(t_start),NGE(t_start),'.','Color',c(dm,:),'MarkerSize',15)
        hold on;
        
    else
        t_start=round(res_lif.MFE_time(i,1)*10);
        t_end=round(res_lif.MFE_time(i+1,1)*10);
        
        plot3(HE(t_start),HI(t_start),NGE(t_start),'.','Color','b','MarkerSize',15)
        hold on;
        
    end
end
% zlim([0 300])
end
