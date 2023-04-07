%% Setting path

addpath('module');

%% parameters for Model_ODE

param.ne       = 300;
param.ni       = 100;
param.s_ee     = 5*0.15;
param.s_ie     = 2*0.5;
param.s_ii     = 4.91*0.4;
param.p_ee     = 1;
param.p_ie     = 1;
param.p_ei     = 1;
param.p_ii     = 1;
param.M        = 100;
param.Mr       = 66;
param.lambda_e = 7;
param.lambda_i = 7;
param.n_exe = 10;
param.tau_ee   = 1.4;
param.tau_ie   = 1.2;
param.tau_ei   = 4.5;
param.tau_ii   = 4.5;
param.duration = 10;
param.delta_time = 0.1;
param.dt = 0.01;
param.bar=50;

%%fig7 plot
figure;
for ei=1:4
subplot(1,4,ei)
param.s_ei     = 1.94+ei*0.03;
init.mean_e=rand(1);
init.var_e=rand(1);
init.mean_i=rand(1);
init.var_i=rand(1);
tic;
res_ode3=model_ode2(param,[]);
toc;
[pos,dm]=detect_MFE(res_ode3);
c=turbo;
ind = length(dm);
HE=res_ode3.h(:,1)*300+res_ode3.h(:,2)*100;
HI=res_ode3.h(:,3)*300+res_ode3.h(:,4)*100;

for i=1:ind-1
    t_start=pos(i);
    t_end=pos(i+1);
    if max(HE(t_start:t_end))<1.5*10^5
    if sum(res_ode3.spikecount_e(t_start:t_end))>100 
        i
        m=ceil(dm(i)*2.3+150);
        plot3(HE(t_start:t_end),HI(t_start:t_end),res_ode3.nge(t_start:t_end),'Color',c(m,:));
         
        hold on;
    else
        plot3(HE(t_start:t_end),HI(t_start:t_end),res_ode3.nge(t_start:t_end),'b');
        
        hold on;
    end
    end
end
view([-30,40])
box on
grid on

xticks([0 60000 120000]);
xticklabels([0 600 1200]);
yticks([0 30000 60000 90000]);
yticklabels([0 300 600 900]);
xlim([0 120000]);
ylim([0 90000]);
title(['S^{EI}=',num2str(1.94+ei*0.03),'×','10^{-2}']);
xlabel('g^{E}')
    ylabel('g^{I}')
    if ei ==1
        zlabel('# of gate E neurons')
    end
end
set(gcf,'Position',[10,10,1600,300]);

%%fig8 plot
dm_all=zeros(50,500);

parfor i=1:50
    i
    param2=param;
    param2.s_ei     = 4.91*(0.4+i*0.001);
    res_ode3=model_ode2(param2,[]);
    [pos,dm]=detect_MFE(res_ode3);
    dm=dm(ceil(length(dm)*0.3):end);
    dm=[length(dm)-1,dm(2:end)-dm(1:end-1),zeros(1,500-length(dm))];
 %   dm=[length(dm),dm(1:end),zeros(1,500-length(dm)-1)];
    dm_all(i,:)=dm;
end

dens=gauss_density(dm_all,[60:-1:-60]);
imagesc(dens);
colormap turbo;
xticks([1 10 20 30 40 50]);
xticklabels([4.91*0.4 4.91*0.41 4.91*0.42 4.91*0.43 4.91*0.44 4.91*0.45])
yticks([1 21 41 61 81 101 121]);
yticklabels([60 40 20 0 -20 -40 -60])
xlabel('S_{ei}');
ylabel('\Deltadm');
set(gcf,'Position',[10,10,1200,400]);
set(gca,'FontSize',16);
set(gca,'FontName','Arial');


%%

function [pos,dm]=detect_MFE(res)
pos=[];
dm=[];
se=res.spikecount_e;
si=res.spikecount_i;
pe=res.peak_e;
pi=res.peak_i;
i=0;
while i<length(se)
    i=i+1;
    if se(i)>0 || (si(i)>0&&sum(si(i-20:i-1))==0)
        j=i-5;
        dm=[dm,p_mean(pe(j,:))-p_mean(pi(j,:))];
        pos=[pos,j];
        i=i+150;
    end
end
end
