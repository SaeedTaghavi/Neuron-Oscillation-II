function []= compare_video(param, res_ode,res_mif, save_path)
hhh=figure;
npis = res_ode.npi;
npes = res_ode.npe;
ve = res_mif.VE;
vi = res_mif.VI;

peak_es = res_ode.peak_e;
peak_is = res_ode.peak_i;
set(gcf,'Position',[10,10,2400,1200]);
ax1 = subplot(2,2,1);
ax2 = subplot(2,2,3);
ax3 = subplot(2,2,2);
ax4 = subplot(2,2,4);
writeObj = VideoWriter([save_path,'ode.avi']);
open(writeObj);
xes = -66:0.01:100;
yes = zeros(1,size(xes,2));
xis = xes;
yis = yes;
for i=1:size(ve,1)
    yes = yes*0;
    npe = npes(i);
    npi = npis(i);
    peak_e = peak_es(i,:);
    peak_i = peak_is(i,:);
    peak_e = reshape(peak_e,[3,10]);
    peak_i = reshape(peak_i,[3,10]);
    x_e = phinv(peak_e(2,1), peak_e(3,1));
    x_e_l = phinv(peak_e(2,1), 0.0013);
    index = (xes>(x_e_l+peak_e(1,1))) &(xes<(x_e + peak_e(1,1)));
    yes(index) = exp(-(xes(index)-peak_e(1,1)).^2./(2*peak_e(2,1)))./(sqrt(2*pi*peak_e(2,1)));
    if npe > 1
        for j=2:npe
            m = peak_e(1,j);
            v = peak_e(2,j);
            r = peak_e(3,j);
            index = (xes>(m-3*sqrt(v))) &(xes<(m+3*sqrt(v)));
            yes(index) = yes(index) + r * exp(-(xes(index)-m).^2./(2*v))./(sqrt(2*pi*v));
        end
    end
    r2 = peak_e(3,npe+1);
    if r2>0
        x_e_u = phinv(peak_e(2,1),peak_e(3,1)+r2);
        index = (xes>0)&(xes<(x_e_u-x_e));
        yes(index) = yes(index) + exp(-(xes(index)+x_e).^2./(2*peak_e(2,1)))./(sqrt(2*pi*peak_e(2,1)));
        
    end
    plot(ax1,xes,yes,'Color','r');
    xlim(ax1,[-66,100]);
    ylim(ax1,[0,0.2]);
    
    yis = yis*0;
    x_i = phinv(peak_i(2,1), peak_i(3,1));
    x_i_l = phinv(peak_i(2,1), 0.0013);
    index = (xis>(x_i_l+peak_i(1,1))) &(xis<(x_i + peak_i(1,1)));
    yis(index) = exp(-(xis(index)-peak_i(1,1)).^2./(2*peak_i(2,1)))./(sqrt(2*pi*peak_i(2,1)));
    if npi > 1
        for j=2:npi
            m = peak_i(1,j);
            v = peak_i(2,j);
            r = peak_i(3,j);
            index = (xis>(m-3*sqrt(v))) &(xis<(m+3*sqrt(v)));
            yis(index) = yis(index) + r * exp(-(xis(index)-m).^2./(2*v))./(sqrt(2*pi*v));
        end
    end
    r2 = peak_i(3,npi+1);
    if r2>0
        x_i_u = phinv(peak_i(2,1),peak_i(3,1)+r2);
        index = (xis>0)&(xis<(x_i_u-x_i));
        yis(index) = yis(index) + exp(-(xis(index)+x_i).^2./(2*peak_i(2,1)))./(sqrt(2*pi*peak_i(2,1)));
    end
    plot(ax2,xis,yis,'Color','b');
    xlim(ax2,[-66,100]);
    ylim(ax2,[0,0.2]);
    
    histogram(ax3,ve(i,:),[-66:1:100],'Normalization','probability');
    xlim(ax3,[-66,100]);
    ylim(ax3,[0,0.2]);
    
    histogram(ax4,vi(i,:),[-66:1:100],'Normalization','probability');
    xlim(ax4,[-66,100]);
    ylim(ax4,[0,0.2]);
    
    sgtitle(['t = ', num2str(i*param.delta_time),' ms']);
    xlabel(ax1,'V_E');
    xlabel(ax2, 'V_I');
    ylabel(ax1,'Percentage');
    ylabel(ax2,'Percentage');
    frame = getframe(hhh);
    writeVideo(writeObj,frame.cdata);
    hold(ax2,'off');
    hold(ax1,'off');
end
close(writeObj);

    function x = phinv(v,h)
        % note: assume original distribution has zero mean so that there is
        % no m variable.
        h=min(h,0.9999);
        h=max(h,0.0001);
        x = erfinv((2*h-1))*sqrt(2)*sqrt(v);
    end
end

