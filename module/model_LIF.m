function [res]=model_LIF(param,init)

%This is the ode version of full model. 
ne = param.ne;
ni = param.ni;
p_ee = param.p_ee;
p_ie = param.p_ie;
p_ei = param.p_ei;
p_ii = param.p_ii;
see = param.s_ee;
sie = param.s_ie;
sei = param.s_ei;
sii = param.s_ii;
s_exe = param.s_exe;
s_exi = param.s_exi;
tau_ee = param.tau_ee;
tau_ie = param.tau_ie;
tau_ei = param.tau_ei;
tau_ii = param.tau_ii;
dt = param.gridsize;

ref_e = param.tau_re;
ref_i = param.tau_ri;

M = param.M;
Mr = param.Mr;
duration = param.duration*1000;

lambda_e = param.lambda_e/1000;
lambda_i = param.lambda_i/1000;
delay = 2;
flag_off_time = 0.0;
start_threshold = 2;
end_threshold = 2;
MFE_interval = 1;
MFE_spike_count_threshold = 5;
flag_time =0;
esi_e=exprnd(1/lambda_e,ceil(lambda_e*(duration*1.4)),ne); %external spike interval
esi_i=exprnd(1/lambda_i,ceil(lambda_i*(duration*1.4)),ni); %external spike interval
while min(sum(esi_e))<duration || min(sum(esi_i))<duration
    esi_e=exprnd(1/lambda_e,ceil(lambda_e*(duration*1.4)),ne); %external spike interval
    esi_i=exprnd(1/lambda_i,ceil(lambda_i*(duration*1.4)),ni); %external spike interval
end



ex_e=zeros(duration/dt,ne);
ex_i=zeros(duration/dt,ni);
%convert esi to the external effect at each time point.
qe=exp(-dt/tau_ee);
qi=exp(-dt/tau_ie);
for i=1:ne
    t=esi_e(1,i);
    count=1;
    while t<duration
        ind=ceil(t/dt);
        ex_e(ind,i)=ex_e(ind,i)+exp((t-ind*dt)/tau_ee)/tau_ee;
        count=count+1;
        t=t+esi_e(count,i);
    end
end
for i=1:ni
    t=esi_i(1,i);
    count=1;
    while t<duration
        ind=ceil(t/dt);
        ex_i(ind,i)=ex_i(ind,i)+exp((t-ind*dt)/tau_ie)/tau_ie;
        count=count+1;
        t=t+esi_i(count,i);
    end
end
for i=2:size(ex_e,1)
    ex_e(i,:)=ex_e(i,:)+ex_e(i-1,:)*qe;
    ex_i(i,:)=ex_i(i,:)+ex_i(i-1,:)*qi;
end

if isempty(init)
    ve=zeros(1,ne);
    vi=zeros(1,ni);
    he=zeros(2,ne);
    hi=zeros(2,ni);
else
    ve=init.ve;
    vi=init.vi;
    he=init.he;
    hi=init.hi;
end

refc_e=zeros(1,ne);%ref clock
refc_i=zeros(1,ni);

sl=duration;
spike_e=zeros(100000,ne);
spike_i=zeros(100000,ni);

nearray=[0:ne-1];
niarray=[0:ni-1];

res.HE = zeros(ceil(duration/dt)+1,ne+ni);
res.HI = zeros(ceil(duration/dt)+1,ne+ni);
res.VE = zeros(ceil(duration/dt)+1,ne);
res.VI = zeros(ceil(duration/dt)+1,ni);
res.spikecount_e = zeros(ceil(duration/dt)+1,1);
res.spikecount_i = zeros(ceil(duration/dt)+1,1);
res.MFE_time = zeros(ceil(duration/dt)+1,2);
res.HEE_stat = zeros(ceil(duration/dt)+1,3);
MFE_time2 = zeros(ceil(duration/dt)+1,2);
HEE_stat2 = zeros(ceil(duration/dt)+1,3);
wave_spike_count2 = zeros(ceil(duration/dt)+1,1);
flag = 0;
wave_spike_count = zeros(ceil(duration/dt)+1,1);
wave_record = zeros(1,100000);
wave_count = 1;
for step=2:duration/dt
    time = step*dt;
    rind_e= refc_e==0; %ref index
    rind_i= refc_i==0;
    refc_e(~rind_e)=refc_e(~rind_e)-1;
    refc_i(~rind_i)=refc_i(~rind_i)-1;
    
    ve(rind_e)=ve(rind_e)+(ex_e(step,rind_e)*s_exe+see*he(1,rind_e)/tau_ee-sei*he(2,rind_e)/tau_ei.*(ve(rind_e)+Mr)/(M+Mr))*dt;
    is_spike=find(0.4*ex_e(step,rind_e)*s_exe<see*he(1,rind_e)/tau_ee & ve(rind_e)>M);
    vi(rind_i)=vi(rind_i)+(ex_i(step,rind_i)*s_exi+sie*hi(1,rind_i)/tau_ie-sii*hi(2,rind_i)/tau_ii.*(vi(rind_i)+Mr)/(M+Mr))*dt;
    
    he=he.*[exp(-dt/tau_ee);exp(-dt/tau_ei)];
    hi=hi.*[exp(-dt/tau_ie);exp(-dt/tau_ii)];
    
    sind_e = ve>M; %spike index
    sind_i = vi>M; %spike index
    spikecount_e=sum(sind_e);
    spikecount_i=sum(sind_i);
    if flag==1
        wave_spike_count(wave_count) = wave_spike_count(wave_count) + spikecount_e + spikecount_i;
    end
    wave_record(wave_record(1)+2: wave_record(1)+ length(is_spike)+1) = time;
    wave_record(1) = wave_record(1) +length(is_spike);
    
    if spikecount_e>0
        if ref_e<0
            ve(sind_e)=ve(sind_e)-M;
        else
            ve(sind_e)=0;
            refc_e(sind_e)=refc_e(sind_e)+ref_e;
        end
        he(1,:)=he(1,:)+binornd(spikecount_e,p_ee,1,ne);
        hi(1,:)=hi(1,:)+binornd(spikecount_e,p_ie,1,ni);
        spike_e(1,sind_e)=spike_e(1,sind_e)+1; %write down the spike time
        spikeind=nearray(sind_e)*100000+spike_e(1,sind_e)+1;
        spike_e(round(spikeind))=step*dt;
    end
    if spikecount_i>0
         if ref_i<0
            vi(sind_i)=vi(sind_i)-M;
         else
            vi(sind_i)=0;
            refc_i(sind_i)=refc_i(sind_i)+ref_i;
         end
        he(2,:)=he(2,:)+binornd(spikecount_i,p_ei,1,ne);
        hi(2,:)=hi(2,:)+binornd(spikecount_i,p_ii,1,ni);
        spike_i(1,sind_i)=spike_i(1,sind_i)+1; %write down the spike time
        spikeind=niarray(sind_i)*100000+spike_i(1,sind_i)+1;
        spike_i(round(spikeind))=step*dt;
    end
    if wave_record(1)>0 && time - wave_record(2) > delay
        num_pop = sum(time - wave_record(2:wave_record(1)+1) >delay);
        wave_record(2:wave_record(1)+1-num_pop) = wave_record(num_pop+2:wave_record(1)+1);
        wave_record(1) = wave_record(1) - num_pop; 
    end
    if flag ==0 && wave_record(1)> start_threshold && time - flag_time >= flag_off_time
        flag = 1;
        flag_time = time;
        res.MFE_time(wave_count,1) = time;
        res.HEE_stat(wave_count,1) = sum(he(1,:));
        HEE_max = sum(he(1,:));
    end
    if flag ==1 && HEE_max < sum(he(1,:))
        HEE_max =sum(he(1,:));
    end
%     if flag ==1 && (wave_record(1)<= end_threshold || (time > flag_time+10 && sum(he(1,:))<50))
    if flag ==1 && (wave_record(1)<= end_threshold || time > flag_time+10 )
        flag = 0;
        flag_time = time;
        res.MFE_time(wave_count, 2) = time;
        res.HEE_stat(wave_count, 2) = HEE_max;
        res.HEE_stat(wave_count, 3) = sum(he(1,:));
        HEE_max = 0;
        wave_count = wave_count +1;
    end
    res.VE(step,:)=ve(:);
    res.VI(step,:)=vi(:);
    res.HE(step,:)=[he(1,:),hi(1,:)];
    res.HI(step,:)=[he(2,:),hi(2,:)];
    res.spikecount_e(step) = spikecount_e;
    res.spikecount_i(step) = spikecount_i;
end
res.spike=[spike_e,spike_i];
res.spike(2:end,:)=res.spike(2:end,:)/1000;

index2 = 1;
index = 1;
while index < wave_count
    MFE_time2(index2,1) = res.MFE_time(index,1);
    HEE_stat2(index2,1) = res.HEE_stat(index,1);
    local_sp_count = 0;
    local_HEE_max = res.HEE_stat(index,2);
    while true
        local_sp_count = local_sp_count + wave_spike_count(index);
        if local_HEE_max < res.HEE_stat(index,2)
            local_HEE_max = res.HEE_stat(index,2);
        end
        HEE_stat2(index2, 3) = res.HEE_stat(index,3);
        MFE_time2(index2, 2) = res.MFE_time(index,2);
        index = index+1;
        if (res.MFE_time(index,1) - res.MFE_time(index-1,2)>= MFE_interval)|| index>wave_count
            break;
        end
    end
    wave_spike_count2(index2) = local_sp_count;
    HEE_stat2(index2,2) = local_HEE_max;
    if local_sp_count >25
        index2  = index2 +1;
    end
end
wave_count = index2;
res.MFE_time = MFE_time2;
res.HEE_stat = HEE_stat2;
wave_spike_count = wave_spike_count2;
res.wave_spike_count = wave_spike_count;
res.wave_count = wave_count;
max_spike = max(res.spike(1,:));
res.spike = res.spike(1:max_spike+1,:);
res.wave_spike_count = res.wave_spike_count(1:wave_count,1);
end


