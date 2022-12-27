%% Setting paths
addpath('module');
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


%% Period 1
param1 = param;
param1.duration = 30;
param1.s_ei = 4.91*0.405/alpha*beta;
tic;
res_lif1=model_LIF3(param1,[]);
toc;

%% Period 2
param2 = param;
param2.duration = 30;
param2.s_ei = 4.91*0.44/alpha*beta;
tic;
res_lif2=model_LIF3(param2,[]);
toc;

%% Period 3
beta = 3;

param3 = param;
param3.duration = 30;
param3.M        = 100*beta;
param3.Mr       = 66*beta;
param3.s_ee     = 5*0.15/alpha*beta;
param3.s_ie     = 2*0.5/alpha*beta;
param3.s_ei     = 4.91*0.418/alpha*beta;
param3.s_ii     = 4.91*0.40/alpha*beta;

param3.lambda_e = 7000*beta;
param3.lambda_i = 7000*beta;
tic;
res_lif3=model_LIF3(param3,[]);
toc;


%% Generate video Beat 1
res = res_lif1;
VE = res.VE/100;
VI = res.VI/100;
timewindow = 200;
num_rec = size(VE,1);
VE = VE(num_rec-timewindow*10: num_rec-1, :);
VI = VI(num_rec-timewindow*10: num_rec-1, :);

vi = VideoWriter('beat-1');
vi.FrameRate = 20;
open(vi);
for i =1:2000
  fi=figure('visible','off');
  subplot(2,1,1);
    a = histogram(VE(i,:));
    a.FaceColor = 'red';
    a.BinWidth = 0.05;
    xlim([-0.67,1]);
    ylim([0,300]);
    ylabel('Number of E neurons');
    subplot(2,1,2);
    b = histogram(VI(i,:));
    b.BinWidth = 0.05;
    b.FaceColor = 'blue';
    xlim([-0.67,1]);
    xlabel('Membrane potential');
    ylim([0,100]);
    ylabel('Number of I neurons');
  sgtitle(['beat 1 t=' num2str(i*0.1, '%0.1f'), ' ms']);
  F = getframe(gcf);
  writeVideo(vi,F);
  close(fi);
end

close(vi);

res = res_lif2;
VE = res.VE/100;
VI = res.VI/100;
timewindow = 200;
num_rec = size(VE,1);
VE = VE(num_rec-timewindow*10: num_rec-1, :);
VI = VI(num_rec-timewindow*10: num_rec-1, :);

vi = VideoWriter('beat-2');
vi.FrameRate = 20;
open(vi);
for i =1:2000
  fi=figure('visible','off');
  subplot(2,1,1);
    a = histogram(VE(i,:));
    a.FaceColor = 'red';
    a.BinWidth =  0.05;
    xlim([-0.67,1]);
    ylim([0,300]);
    ylabel('Number of E neurons');
    subplot(2,1,2);
    b = histogram(VI(i,:));
    b.BinWidth = 0.05;
    b.FaceColor = 'blue';
    xlim([-0.67,1]);
    xlabel('Membrane potential');
    ylim([0,100]);
    ylabel('Number of I neurons');
  sgtitle(['beat 2 t=', num2str(i*0.1, '%0.1f'), ' ms']);
  F = getframe(gcf);
  writeVideo(vi,F);
  close(fi);
end

close(vi);

res = res_lif3;
VE = res.VE/300;
VI = res.VI/300;
timewindow = 200;
num_rec = size(VE,1);
VE = VE(num_rec-timewindow*10: num_rec-1, :);
VI = VI(num_rec-timewindow*10: num_rec-1, :);

vi = VideoWriter('beat-3');
vi.FrameRate = 20;
open(vi);
for i =1:2000
  fi=figure('visible','off');
  subplot(2,1,1);
    a = histogram(VE(i,:));
    a.FaceColor = 'red';
    a.BinWidth =  0.05;
    xlim([-0.67,1]);
    ylim([0,300]);
    ylabel('Number of E neurons');
    subplot(2,1,2);
    b = histogram(VI(i,:));
    b.BinWidth = 0.05;
    b.FaceColor = 'blue';
    xlim([-0.67,1]);
    xlabel('Membrane potential');
    ylim([0,100]);
    ylabel('Number of I neurons');
  sgtitle(['beat 3 t=' num2str(i*0.1, '%0.1f'), ' ms']);
  F = getframe(gcf);
  writeVideo(vi,F);
  close(fi);
end

close(vi);


