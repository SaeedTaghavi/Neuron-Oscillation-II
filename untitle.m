p = 0.8;
param.s_exe = 1;
param.s_exi = 1;
param2 = param;
param2.p_ee     = p;
param2.p_ie     = p;
param2.p_ei     = p;
param2.p_ii     = p;
param2.s_ee     = 5*0.15/p;
param2.s_ie     = 2*0.5/p;
param2.s_ei     = 4.91*0.5/p;
param2.s_ii     = 4.91*0.4/p;

ve = zeros(1,300);
vi = zeros(1,100);
he = zeros(2,300);
hi = zeros(2,100);

%%
res_lif = model_LIF3(param2, ve, vi, he, hi);

%%
figure;
rasterplot(res_lif, param2);
hold on;
for i =1:res_lif.wave_count-1
    xline(res_lif.MFE_time(i,1),'b','LineWidth',1);
    hold on;
    xline(res_lif.MFE_time(i,2),'LineWidth',1);
    hold on;
end
