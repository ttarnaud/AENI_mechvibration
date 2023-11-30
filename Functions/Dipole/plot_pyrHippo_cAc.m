% import data
clear all
close all
clc
Gif_flag = 0;
sim_dip = 7;
Plot_V_flag = 0;
Interp_flag =  [1, 1];
switch sim_dip
    case 1
        filename = 'dipole_pyr_cAC';
        dataLines = [4, 5491];        
        dipole_data.V = import_V_dipole_neuron(filename,dataLines);
        dataLines2 = [5496+261, 11242];
        Plot_V_flag = 1;
        Fs1 = 40e3;
        
    case 2
        filename = 'hippocampus_dipole_pyr_cAC_2_fixedstep-5us';
        dataLines2 = [4, 20005];
        Fs1 = 400e3;
        Interp_flag(1) = 0;
        Majortitle = 'pyramidal hippocampal neuron (2), step = 5탎';
    case 3
        filename = 'hippocampus_dipole_pyr_cAC_2_fixedstep-25us';
        dataLines2 = [4, 4004];
        Fs1 = 40e3;
        Interp_flag(1) = 0;
        Majortitle = 'pyramidal hippocampal neuron (2), step = 25탎';
    case 4
        filename = 'hippocampus_dipole_pyr_cAC_fixedstep-5us';
        dataLines2 = [4, 20005];
        Fs1 = 400e3;
        Interp_flag(1) = 0;
        Majortitle = 'pyramidal hippocampal neuron, step = 5탎';
    case 5
        filename = 'hippocampus_dipole_pyr_cAC_fixedstep-25us';
        dataLines2 = [4, 4004];
        Fs1 = 40e3;
        Interp_flag(1) = 0;
        Majortitle = 'pyramidal hippocampal neuron, step = 25탎';
    case 6
        filename = 'L23_PC_cADpyr229_3_dipole_fixedstep-5us';
        dataLines2 = [4, 20003];
        Fs1 = 400e3;
        Interp_flag(1) = 0;
        Majortitle = 'pyramidal cortex L23 neuron, step = 5탎';
    case 7
        filename = 'L23_PC_cADpyr229_3_dipole_fixedstep-25us';
        dataLines2 = [4, 4003];
        Fs1 = 40e3;
        Interp_flag(1) = 0;
        Majortitle = 'pyramidal cortex L23 neuron, step = 25탎';
        
end

fmiddle = 35e3;

dipole_data.Min = import_dipole_moment(filename,dataLines2);

[dipole_data.M.Time,ia] = unique([dipole_data.Min.Time]);
dipole_data.M.dipx = dipole_data.Min.dipx(ia);
dipole_data.M.dipy = dipole_data.Min.dipy(ia);
dipole_data.M.dipz = dipole_data.Min.dipz(ia);
dipole_data.M.diptot = dipole_data.Min.diptot(ia);

%% Plot
figure
if Plot_V_flag
subplot(2,1,1)
plot(dipole_data.V.Time,dipole_data.V.V_soma)
xlabel('Time [ms]')
ylabel('V_{soma} [mV]')
subplot(2,1,2)
end

plot(dipole_data.M.Time, dipole_data.M.dipx,'DisplayName','dip_x')
hold on
plot(dipole_data.M.Time, dipole_data.M.dipy,'DisplayName','dip_y')
plot(dipole_data.M.Time, dipole_data.M.dipz,'DisplayName','dip_z')
plot(dipole_data.M.Time, dipole_data.M.diptot,'DisplayName','dip_{tot}')
hold off
%% rotate dominant axis onto y
x_avg = mean(abs(dipole_data.M.dipx));
y_avg = mean(abs(dipole_data.M.dipy));
z_avg = mean(abs(dipole_data.M.dipz));
xline = [-10:0.01:10]*x_avg;
yline = [-10:0.01:10]*y_avg;
zline = [-10:0.01:10]*z_avg;

yaxis = [0,1,0];
dominant = [x_avg,y_avg,z_avg]; dominant =  dominant/norm(dominant);

v = cross(dominant,yaxis);
s = norm(v); %sin of angle
c = dot(yaxis,dominant); %cosine of angle
skewv = [0 -v(3) v(2);v(3) 0 -v(1);-v(2) v(1) 0];
R = eye(3)+skewv+skewv^2.*(1-c)/s^2;

newdipole = R*[dipole_data.M.dipx';dipole_data.M.dipy';dipole_data.M.dipz'];
newdipolex = newdipole(1,:);
newdipoley = newdipole(2,:);
newdipolez = newdipole(3,:);
newdipoletot = vecnorm(newdipole,2,1);
%% create periodograms

Fss = [Fs1,40e6];
figure(10)
for iplot = 1:2
    Fs = Fss(iplot);
    if Interp_flag(iplot)
    
t1 = dipole_data.M.Time(1)*1e-3:1/Fs:dipole_data.M.Time(end)*1e-3;
dipx_inter = interp1(dipole_data.M.Time*1e-3,dipole_data.M.dipx,t1);
dipy_inter = interp1(dipole_data.M.Time*1e-3,dipole_data.M.dipy,t1);
dipz_inter = interp1(dipole_data.M.Time*1e-3,dipole_data.M.dipz,t1);
diptot_inter = interp1(dipole_data.M.Time*1e-3,dipole_data.M.diptot,t1);

newdipx_inter = interp1(dipole_data.M.Time*1e-3,newdipolex,t1);
newdipy_inter = interp1(dipole_data.M.Time*1e-3,newdipoley,t1);
newdipz_inter = interp1(dipole_data.M.Time*1e-3,newdipolez,t1);
newdiptot_inter = interp1(dipole_data.M.Time*1e-3,newdipoletot,t1);
    else
        dipx_inter = dipole_data.M.dipx;
        dipy_inter = dipole_data.M.dipy;
        dipz_inter = dipole_data.M.dipz;
        diptot_inter = dipole_data.M.diptot;
        
        newdipx_inter = newdipolex;
        newdipy_inter = newdipoley;
        newdipz_inter = newdipolez;
        newdiptot_inter = newdipoletot;
    end

[PSDdipx,fx] = periodogram(dipx_inter,rectwin(length(dipx_inter)),length(dipx_inter),Fs,'power');
[PSDdipy,fy] = periodogram(dipy_inter,rectwin(length(dipy_inter)),length(dipy_inter),Fs,'power');
[PSDdipz,fz] = periodogram(dipz_inter,rectwin(length(dipz_inter)),length(dipz_inter),Fs,'power');
[PSDdiptot,ftot] = periodogram(diptot_inter,rectwin(length(diptot_inter)),length(diptot_inter),Fs,'power');

[nPSDdipx,fx] = periodogram(newdipx_inter,rectwin(length(newdipx_inter)),length(newdipx_inter),Fs,'power');
[nPSDdipy,fy] = periodogram(newdipy_inter,rectwin(length(newdipy_inter)),length(newdipy_inter),Fs,'power');
[nPSDdipz,fz] = periodogram(newdipz_inter,rectwin(length(newdipz_inter)),length(newdipz_inter),Fs,'power');
[nPSDdiptot,ftot] = periodogram(newdiptot_inter,rectwin(length(newdiptot_inter)),length(newdiptot_inter),Fs,'power');

figure(10)
subplot(2,2,2*(iplot-1)+1)
plot(ftot,10*log10(PSDdiptot),'DisplayName','diptot');
hold on
plot(fx,10*log10(PSDdipx),'DisplayName','dipx')
plot(fy,10*log10(PSDdipy),'DisplayName','dipy');
plot(fz,10*log10(PSDdipz),'DisplayName','dipz');

hold off
set(gca,'xscale','log')%,'yscale','log')
ylabel('Power [dB]')
if iplot == 1
    title('original')
end
if iplot==2
    xlabel('Frequency [Hz]')
end

subplot(2,2,2*(iplot-1)+2)
plot(ftot,10*log10(nPSDdiptot),'DisplayName','diptot')
hold on
plot(fy,10*log10(nPSDdipy),'DisplayName','dipy')
plot(fx,10*log10(nPSDdipx),'DisplayName','dipx')
plot(fz,10*log10(nPSDdipz),'DisplayName','dipz')

legend('show')
hold off
set(gca,'xscale','log')%,'yscale','log')
if iplot == 1
    title('transformed y=principal component')
end
if iplot==2
    xlabel('Frequency [Hz]')
end

n_ndpy = nPSDdipy./nPSDdipy(1);
% isolate first 3 decades
figure(20)
subplot(1,2,iplot)
p1 = plot(fy,10*log10(n_ndpy));
hold on
set(get(get(p1,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');

f_idxlow = 10<fy & fy<1e3;
n_ndpy_tokHz = n_ndpy(f_idxlow);

ab_low = [ones(size(fy(f_idxlow))),log10(fy(f_idxlow))]\(10*log10(n_ndpy_tokHz));
plot(fy(f_idxlow),ab_low(1)+log10(fy(f_idxlow))*ab_low(2),'DisplayName',['alpha = ',num2str(round(ab_low(2)/10,1))],...
    'LineWidth',2)

% isolate rest to sampling frequency
f_idxmiddle = 1e3<fy & fy<=fmiddle;
n_ndpy_toFsmiddle = n_ndpy(f_idxmiddle);
ab_middle = [ones(size(fy(f_idxmiddle))),log10(fy(f_idxmiddle))]\(10*log10(n_ndpy_toFsmiddle));
plot(fy(f_idxmiddle),ab_middle(1)+log10(fy(f_idxmiddle))*ab_middle(2),'DisplayName',['alpha = ',num2str(round(ab_middle(2)/10,1))],...
    'LineWidth',2)

% isolate rest to sampling frequency
f_idxhigh = fmiddle<fy & fy<=Fss(2);
n_ndpy_toFs = n_ndpy(f_idxhigh);
ab_high = [ones(size(fy(f_idxhigh))),log10(fy(f_idxhigh))]\(10*log10(n_ndpy_toFs));
plot(fy(f_idxhigh),ab_high(1)+log10(fy(f_idxhigh))*ab_high(2),'DisplayName',['alpha = ',num2str(round(ab_high(2)/10,1))],...
    'LineWidth',2)
hold off
legend('show')
set(gca,'xscale','log')
if Interp_flag(iplot)
    title_str = ['interpolated to Fs =  ', num2str(Fss(iplot), '%3.2e')];
else
    title_str = ['not interpolated, Fs = ', num2str(Fss(iplot), '%3.2e')];
end
title(title_str)
Ylim = get(gca,'ylim');
set(gca,'ylim',[Ylim(1),0])
if iplot==1
ylabel({Majortitle,'Power [dB]'})
end
xlabel('Frequency [Hz]')
set(findall(gcf,'-property','FontSize'),'FontSize',14)

end
set(gcf,'position',[380,557,1324,421]);
figure(10)
set(gcf,'position',[380, 163,1239, 815]);
mtit(Majortitle,'xoff',0,'yoff',.05);

%%

% 
% figure
% plot(fy,10*log10(nPSDdipy))
% hold on
% plot(fy,10*log10(n_ndpy))
% hold off
% set(gca,'xscale','log')




%% create giffs
if Gif_flag
L = length(dipole_data.M.Time);
dipolemoment = figure;
Gifname = 'testAnimatedDipole.gif';
Step = ceil(L/100);
DelayTime = Step*dipole_data.M.Time(end)/2000;

viewsp = [67.5,30];

Xlims = 1.1*[min(dipole_data.M.dipx),max(dipole_data.M.dipx)];
Ylims = 1.1*[min(dipole_data.M.dipy),max(dipole_data.M.dipy)];
Zlims = 1.1*[min(dipole_data.M.dipz),max(dipole_data.M.dipz)];
Ylim2 = 1.1*[-max(dipole_data.M.diptot),max(dipole_data.M.diptot)];

Xlims3 = 1.1*[min(newdipolex),max(newdipolex)];
Ylims3 = 1.1*[min(newdipoley),max(newdipoley)];
Zlims3 = 1.1*[min(newdipolez),max(newdipolez)];


x_avg = mean(abs(dipole_data.M.dipx));
y_avg = mean(abs(dipole_data.M.dipy));
z_avg = mean(abs(dipole_data.M.dipz));
xline = [-10:0.01:10]*x_avg;
yline = [-10:0.01:10]*y_avg;
zline = [-10:0.01:10]*z_avg;

xline2 = [-10:0.01:10]*0;
yline2 = [-10:0.01:10]*1;
zline2 = [-10:0.01:10]*0;
for i=1:Step:L
    figure(dipolemoment)
    set(gcf,'position',[-1695,87, 1350,851])
    X = [0, dipole_data.M.dipx(i)];
    Y = [0, dipole_data.M.dipy(i)];
    Z = [0, dipole_data.M.dipz(i)];
    X2 = X; Y2 = Y; Z2 = Z;
    if X(2)<X(1)
        X2 = abs(fliplr(X));
    end
    if Y(2)<Y(1)
        Y2 = abs(fliplr(Y));
    end
    if Z(2)<Z(1)
        Z2 = abs(fliplr(Z));
    end
    radius = dipole_data.M.diptot(i)/10;
    
    Xnew = [0, newdipolex(i)];
    Ynew = [0, newdipoley(i)];
    Znew = [0, newdipolez(i)];
    Xnew2 = Xnew; Ynew2 = Ynew; Znew2 = Znew;
    if Xnew(2)<Xnew(1)
        Xnew2 = abs(fliplr(Xnew));
    end
    if Ynew(2)<Ynew(1)
        Ynew2 = abs(fliplr(Ynew));
    end
    if Znew(2)<Znew(1)
        Znew2 = abs(fliplr(Znew));
    end
    radiusnew = newdipoletot(i)/10;
    
%     subplot(2,3,1)
%     plot3(nan,nan,nan)
%     hold on
%     xlim(Xlims)
%     ylim(Ylims)
%     zlim(Zlims)
%     xlabel('x-axis')
%     ylabel('y-axis')
%     zlabel('z-axis')
%     title(['Time: ', num2str(dipole_pyr_cAC.M.Time(i))])
%     plot3(xline,yline,zline,'r')
%     arrow3d(X,Y,Z,0.75,radius);
%     hold off
    

    subplot(3,3,1)
    plot3(nan,nan,nan)
    hold on
    xlim(Xlims)
    ylim(Ylims)
    zlim(Zlims)
    xlabel('X-axis')
    ylabel('Y-axis')
    zlabel('Z-axis')
    title(['Time: ', num2str(dipole_data.M.Time(i))])
    %set(gca,'view',viewsp)
    plot3(xline,yline,zline,'r')
    arrow3d(X2,Y2,Z2,0.75,radius);
    hold off
    viewsp2 = get(gca,'view');
    
    % plot Y direction up so permute plot info Y on Z Z on X en X on Y
    subplot(3,3,2)
    hold on
    xlim(Xlims)
    ylim(Ylims)
    zlim(Zlims)
    xlabel('X-axis')
    ylabel('Y-axis')
    zlabel('Z-axis')
    set(gca,'view',viewsp2)
    title(['Time: ', num2str(dipole_data.M.Time(i))])
    plot3(xline,yline,zline,'r')
    arrw = arrow3d(X2,Y2,Z2,0.75,radius);
    set(arrw,'facealpha',0.1)
    hold off
    
    % plot Y direction up so permute plot info Y on Z Z on X en X on Y
    subplot(3,3,3)
    plot3(nan,nan,nan)
    hold on
    xlim(Xlims3)
    ylim(Ylims3)
    zlim(Zlims3)
    xlabel('X-axis')
    ylabel('Y-axis')
    zlabel('Z-axis')
    title(['Projected Time: ', num2str(dipole_data.M.Time(i))])
    set(gca,'view',viewsp)
    
    plot3(xline2,yline2,zline2,'r')
    arrow3d(Xnew2,Ynew2,Znew2,0.75,radiusnew);
    hold off
    
    
    subplot(3,3,[4,5,6])
    plot(nan,nan,'HandleVisibility','off')
    hold on
    xlim([0,dipole_data.M.Time(end)])
    ylim(Ylim2);
    xlabel('Time [ms]')
    plot(dipole_data.M.Time(1:i),dipole_data.M.dipx(1:i),'DisplayName','dip_x')
    plot(dipole_data.M.Time(1:i),dipole_data.M.dipy(1:i),'DisplayName','dip_y')
    plot(dipole_data.M.Time(1:i),dipole_data.M.dipz(1:i),'DisplayName','dip_z')
    plot(dipole_data.M.Time(1:i),dipole_data.M.diptot(1:i),'DisplayName','dip_{tot}')
    legend('show')
    hold off
    title('Output Neuron')
    drawnow
    
    subplot(3,3,[4,5,6]+3)
    plot(nan,nan,'HandleVisibility','off')
    hold on
    xlim([0,dipole_data.M.Time(end)])
    ylim(Ylim2);
    xlabel('Time [ms]')
    plot(dipole_data.M.Time(1:i),newdipolex(1:i),'DisplayName','dip_x')
    plot(dipole_data.M.Time(1:i),newdipoley(1:i),'DisplayName','dip_y')
    plot(dipole_data.M.Time(1:i),newdipolez(1:i),'DisplayName','dip_z')
    plot(dipole_data.M.Time(1:i),newdipoletot(1:i),'DisplayName','dip_{tot}')
    legend('show')
    hold off
    title('After projection on dominant direction')
    drawnow
    
    pause(0.001)
    % Capture the plot as an image 
    frame = getframe(dipolemoment);
    im = frame2im(frame);
    [imind,cm] = rgb2ind(im,256);
    
    % Write to the GIF File 
      if i == 1 
          imwrite(imind,cm,Gifname,'gif', 'Loopcount',inf,'DelayTime',DelayTime); 
      else 
          imwrite(imind,cm,Gifname,'gif','WriteMode','append','DelayTime',DelayTime); 
      end 
    
end
end


