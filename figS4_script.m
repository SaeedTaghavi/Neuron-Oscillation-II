%% figS4
% The distribution of the states of MFE initiations in the reduced state space. 
%% Setting Path

addpath('module');


%% Parameters


taus = 0:0.5:5.5;
sexs = 0.2:0.8:9;
ptauEs = [0.75,1,2,3];
ps=0.25:0.25:1;
seis=2.44:0.08:2.68;
tod = datetime('now','Format','yyMMdd');
tod = char(tod);
p=0.8;
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
param.s_ii     = 4.91*0.40/alpha*beta/p;

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
param.lambda_e = 7000*beta/sexs(2);
param.lambda_i = 7000*beta/sexs(2);
param.gridsize=0.1;

%% Models with Change S_EI

NGEs = cell(12,1);
NGIs = cell(12,1);
HEs = cell(12,1);
HIs = cell(12,1);
RESs = cell(12,1);
MEs = cell(12,1);
MIs = cell(12,1);

for i=1:4
    param2 = param;
    param2.s_ei = seis(i);
    
    tic;
    res_lif=model_LIF_FC(param2,[]);
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
end

%% Subfig with change S_EI
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
    xticks([50 55 60 65 70 75 80 85 90]);
    yticks([30 40 50 60 70 80 90])

    zticks([0 2000 4000 6000]);
    zticklabels([0 20 40 60]);

     title(['S^{EI}=',num2str(seis(i),'%4.2f'),'×','10^{-2}']);
    view([-25,55]);
%     xlim([65 90]);
%     ylim([50 90]);
    zlim([0 6000]);

    set(gca,'fontsize',15,'fontname','Arial'); 
end
set(gcf,'Position',[10,10,1600,300]);
f = gcf;
exportgraphics(f,['figs/FigS4/A-', tod,'.eps'],'Resolution',300)

%% Models with Change S_ext 
NGEs = cell(12,1);
NGIs = cell(12,1);
HEs = cell(12,1);
HIs = cell(12,1);
RESs = cell(12,1);
MEs = cell(12,1);
MIs = cell(12,1);

for i=1:4
    param2 = param;
    param2.s_exe = sexs((i-1)*3+1);
    param2.s_exi = sexs((i-1)*3+1);
    param2.lambda_e = 7000*beta/sexs((i-1)*3+1);
    param2.lambda_i = 7000*beta/sexs((i-1)*3+1);
    % Model LIF

    tic;
    res_lif=model_LIF_FC(param2,[]);
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
end
c=turbo;

%% Subfig with change S_ext
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
    xticks([50 55 60 65 70 75 80 85 90]);
    yticks([30 40 50 60 70 80 90])

    zticks([0 2000 4000 6000]);
    zticklabels([0 20 40 60]);
         title(['S^{ext}=',num2str(sexs((i-1)*3+1)),'×','10^{-2}']);
    view([-25,55]);
    zlim([0 6000]);
%     xlim([50 90]);
%     ylim([50 90]);
    set(gca,'fontsize',15,'fontname','Arial'); 
end
set(gcf,'Position',[10,10,1600,300]);
f = gcf;
exportgraphics(f,['figs/FigS4/B-', tod,'.eps'],'Resolution', 300);

%% Change P 
NGEs = cell(12,1);
NGIs = cell(12,1);
HEs = cell(12,1);
HIs = cell(12,1);
RESs = cell(12,1);
MEs = cell(12,1);
MIs = cell(12,1);

for i=1:4
    param2 = param;
    p=ps(i);
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
end
c=turbo;

%% Fig of Change P 
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
    xticks([50 55 60 65 70 75 80 85 90]);
    yticks([30 40 50 60 70 80 90])

    zticks([0 2000 4000 6000]);
    zticklabels([0 20 40 60]);
    title(['P=',num2str(ps(i),'%4.2f')]);
    view([-25,55]);
    zlim([0 6000]);
%     xlim([50 90]);
%     ylim([50 90]);
    set(gca,'fontsize',15,'fontname','Arial'); 
end
set(gcf,'Position',[10,10,1600,300]);
f = gcf;
exportgraphics(f,['figs/FigS4/C-', tod,'.eps'],'Resolution', 300);

%% Change Tau R
NGEs = cell(12,1);
NGIs = cell(12,1);
HEs = cell(12,1);
HIs = cell(12,1);
RESs = cell(12,1);
MEs = cell(12,1);
MIs = cell(12,1);

for i=1:4
    param2 = param;
    param2.tau_re=taus((i-1)*3+1);
    param2.tau_ri=taus((i-1)*3+1);

    % Model LIF

    tic;
    res_lif=model_LIF_FC(param2,[]);
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
end
c=turbo;

%% Fig of Change tau_R 
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
    xticks([50 55 60 65 70 75 80 85 90]);
    yticks([30 40 50 60 70 80 90]);

    zticks([0 2000 4000 6000]);
    zticklabels([0 20 40 60]);
    title(['\tau_R=',num2str(taus((i-1)*3+1),'%4.1f')]);
    view([-25,55]);
    zlim([0 6000]);
%     xlim([50 90]);
%     ylim([50 90]);
    set(gca,'fontsize',15,'fontname','Arial'); 
end
set(gcf,'Position',[10,10,1600,300]);
f = gcf;
exportgraphics(f,['figs/FigS4/D-', tod,'.eps'],'Resolution', 300);

%% Change tau_E

NGEs = cell(12,1);
NGIs = cell(12,1);
HEs = cell(12,1);
HIs = cell(12,1);
RESs = cell(12,1);
MEs = cell(12,1);
MIs = cell(12,1);

for i=1:4
    param2 = param;
    param2.tau_ee = 1.4*ptauEs(i);
    param2.tau_ie = 1.2*ptauEs(i);

    % Model LIF

    tic;
    res_lif=model_LIF_FC(param2,[]);
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
end
c=turbo;


%% FIG of change Tau_E
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
    xticks([50 55 60 65 70 75 80 85 90]);
    yticks([30 40 50 60 70 80 90]);

    zticks([0 2000 4000 6000 8000 10000]);
    zticklabels([0 20 40 60 80 100]);
title(['\tau_E =', num2str(ptauEs(i)),' times orginal values']);
    view([-25,55]);
    zlim([0 10000]);
%     xlim([50 90]);
%     ylim([50 90]);
    set(gca,'fontsize',15,'fontname','Arial'); 
end
set(gcf,'Position',[10,10,1600,300]);
f = gcf;
exportgraphics(f,['figs/FigS4/E-', tod,'.eps'],'Resolution', 300);

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
 zlim([0 300]);
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
