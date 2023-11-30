% dipolemoment figure in appendix
clear all
close all
clc

sim_dps = 3:2:7;
idx_sp = 0;
fignr = 10;
psd_color = [0.4,0.4,0.4];
psd_low_color = [228,26,28]/256;%[1,0.3,0.3];
psd_middle_color = [55,126,184]/256;%[0.3,1,0.3];
psd_high_color = [77,175,74]/256;%[0.3,0.3,1];
fontsize = 14;
flow = 500;
fmiddle = 5000;%25e3;

l_PSD = {};
for sim_dip = sim_dps

idx_sp =idx_sp+1;


switch sim_dip
    case 1
        filename = 'dipole_pyr_cAC';
        dataLines = [4, 5491];        
        dipole_data.V = import_V_dipole_neuron(filename,dataLines);
        dataLines2 = [5496+261, 11242];
        Fs1 = 40e3;
        
        
    case 2
        filename = 'hippocampus_dipole_pyr_cAC_2_fixedstep-5us';
        dataLines2 = [4, 20005];
        Fs1 = 400e3;
        Majortitle = 'pyramidal hippocampal neuron (2), step = 5탎';
    case 3
        filename = 'hippocampus_dipole_pyr_cAC_2_fixedstep-25us'; %cell morpho mpg141209_A_idA.asc
        dataLines2 = [4, 4004];
        Fs1 = 40e3;
        Majortitle = 'pyramidal hippocampal neuron (2), step = 25탎';
        label = 'hip 1209\_A';
    case 4
        filename = 'hippocampus_dipole_pyr_cAC_fixedstep-5us'; % cell morpho mpg141208_B_idA.asc
        dataLines2 = [4, 20005];
        Fs1 = 400e3;
        Majortitle = 'pyramidal hippocampal neuron, step = 5탎';
    case 5
        filename = 'hippocampus_dipole_pyr_cAC_fixedstep-25us'; %cell morpho mpg141208_B_idA.asc
        dataLines2 = [4, 4004];
        Fs1 = 40e3;
        Majortitle = 'pyramidal hippocampal neuron, step = 25탎';
        label = 'hip 1208\_B';
    case 6
        filename = 'L23_PC_cADpyr229_3_dipole_fixedstep-5us';
        dataLines2 = [4, 20003];
        Fs1 = 400e3;
        Majortitle = 'pyramidal cortex L23 neuron, step = 5탎';
    case 7
        filename = 'L23_PC_cADpyr229_3_dipole_fixedstep-25us';
        dataLines2 = [4, 4003];
        Fs1 = 40e3;
        Majortitle = 'pyramidal cortex L23 neuron, step = 25탎';
        label = 'L23 pyr 229\_3';
        
end
%Load data
dipole_data.Min = import_dipole_moment(filename,dataLines2);

[time,ia] = unique([dipole_data.Min.Time]);
dipx = dipole_data.Min.dipx(ia);
dipy = dipole_data.Min.dipy(ia);
dipz = dipole_data.Min.dipz(ia);
diptot = dipole_data.Min.diptot(ia);

%Calculate PSD
[PSDdipx,fx] = periodogram(dipx,rectwin(length(dipx)),length(dipx),Fs1,'power');
[PSDdipy,fy] = periodogram(dipy,rectwin(length(dipy)),length(dipy),Fs1,'power');
[PSDdipz,fz] = periodogram(dipz,rectwin(length(dipz)),length(dipz),Fs1,'power');
[PSDdiptot,ftot] = periodogram(diptot,rectwin(length(diptot)),length(diptot),Fs1,'power');

%normalize
fi = ftot;
PSDi = PSDdiptot;
nPSD = PSDdiptot/PSDdiptot(1);

%calculate powerlaws
f_idxlow = 0<fi & fi<flow;
nPSD_tokHz = nPSD(f_idxlow);
ab_low = [ones(size(fi(f_idxlow))),log10(fi(f_idxlow))]\(log10(nPSD_tokHz));
PSD_low = 10*(ab_low(1)+log10(fi(f_idxlow))*ab_low(2));

f_idxmiddle = flow<fi & fi<=fmiddle;
mPSD_toFsmiddle = nPSD(f_idxmiddle);
ab_middle = [ones(size(fi(f_idxmiddle))),log10(fi(f_idxmiddle))]\(log10(mPSD_toFsmiddle));
PSD_middle = 10*(ab_middle(1)+log10(fi(f_idxmiddle))*ab_middle(2));

f_idxhigh = fmiddle<fi;
nPSD_toFs = nPSD(f_idxhigh);
ab_high = [ones(size(fi(f_idxhigh))),log10(fi(f_idxhigh))]\(log10(nPSD_toFs));
PSD_high = 10*(ab_high(1)+log10(fi(f_idxhigh))*ab_high(2));


% Plot
figure(fignr)

% time plot
subplot(length(sim_dps),2,(idx_sp-1)*2+1)
plot(time, dipx,'DisplayName','dp_x')
hold on
plot(time, dipy,'DisplayName','dp_y')
plot(time, dipz,'DisplayName','dp_z')
plot(time, diptot,'DisplayName','| dp |')
hold off
if idx_sp == length(sim_dps); xlabel('time [ms]','fontsize',fontsize,'Interpreter','latex'); end
if idx_sp == 1; ti_time = title('dp [pA$\cdot$m]','Interpreter','latex','fontsize',fontsize);end
ylabel(label,'Interpreter','latex')
set(gca,'box','off')
if idx_sp == length(sim_dps); legend_time = legend('show','box','off','NumColumns',4,'location','southoutside'); end

%PSD plot
subplot(length(sim_dps),2,idx_sp*2)
p1 = plot(fi,10*log10(nPSD),'color',psd_color,'HandleVisibility', 'off');
if idx_sp == 1;ti_PSD = title('Power [dB]','Interpreter','latex','fontsize',fontsize);end
hold on
plot(fi(f_idxlow),PSD_low,'color',psd_low_color,'DisplayName',['$\alpha_{l}$ = ',num2str(round(ab_low(2),1))],...
    'LineWidth',2)

plot(fi(f_idxmiddle),PSD_middle,'color',psd_middle_color,'DisplayName',['$\alpha_{m}$ = ',num2str(round(ab_middle(2),1))],...
    'LineWidth',2)

plot(fi(f_idxhigh),PSD_high,'color',psd_high_color,'DisplayName',['$\alpha_{h}$ = ',num2str(round(ab_high(2),1))],...
    'LineWidth',2)
hold off
l_PSD{idx_sp} = legend('show','box','off','Interpreter','latex');
set(gca,{'xscale','box','xtick'},{'log','off',[1e0,1e1,1e2,1e3,1e4,1e5]})


end
xlabel('frequency [Hz]', 'Interpreter','latex')
set(legend_time,'position',[0.0698    0.9011    0.4398    0.0481])
set(ti_time,'position',(get(ti_time,'position').*[1,1.1,1]))
set(ti_PSD,'position',(get(ti_PSD,'position').*[1,2,1]))
for i = 1:length(l_PSD)
    set(l_PSD{i},'position',(get(l_PSD{i},'position')+[0.07,0.032,0,0]));%.*[1.09,1.04,1,1]));
end
set(findall(gcf,'type','axes'),'fontsize',11)
set(gcf,{'units','color','position','paperunits','papersize','Renderer'},...
    {'centimeters',[1,1,1],[1,3,17.8,13.2],'centimeters',[17.8,13.2+1],'Painters'})