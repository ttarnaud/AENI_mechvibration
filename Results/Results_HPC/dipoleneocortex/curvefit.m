%curve fit and plot extension
clear all
close all
addpath(pwd)
%% results with three models
Filenames = dir('Results_*');
Filenames = Filenames(2:end); % to exclude settingset13 old
for ifnms = 1:length(Filenames)
name = Filenames(ifnms).name(9:end-4);
load(['Results_',name,'.mat'])
open(['Strngthcurves',name,'.fig'])

h = findobj(gcf,'type','axes');

extrapol = cell(1,4);
extrapol{1} = name(1:end-14);
for imodel = 1:3
ydata = squeeze(requiredAmp(imodel,:,:))*1e9;
xdata = logspace(0,floor(length(ydata)/2),length(ydata))';
[fitresult,coeff,gof] = FitPower(xdata,ydata);
coeff = round(coeff,2,'significant');
fun = @(x) coeff(1).*x.^coeff(2)+coeff(3);
namefun = ['fit: ',num2str(coeff(1)),' x^{',num2str(coeff(2)),'} + ',num2str(coeff(3))];
axes(h(4-imodel))
set(findobj(gca,'type','line'),'LineStyle','none')
set(findobj(gca,'type','line'),'Marker','*')
hold on
plot(logspace(0,5,101),fun(logspace(0,floor(length(ydata)/2),101)),'DisplayName',namefun)
hold off
if imodel == 1
    legend({'RMS threshold = 10%',namefun})
end
if imodel == 3
    extrapol{1+imodel} = fun(10^4);
else
    extrapol{1+imodel} = fun(10^6);
end
end
if ifnms == 1
% Convert cell to a table and use first row as variable names
T = cell2table(extrapol,'VariableNames',{'SettingSet','HumansScalp','HumanCortical','MouseScalp'});
else
    T = [T;extrapol];
end

end
%% results of single runs
OLDdirectory = cd('./Singleruns');
Filenames = dir('Results_*');
missedSettingsets = {'Settingset3','Settingset9','Settingset10','Settingset11','Settingset12'};
Titles = {'HumanScalp','HumanCortical','MouseScalp'};
for ifnms = 1:length(Filenames)
name = Filenames(ifnms).name(9:end-4);
load(['Results_',name,'.mat'])
if mod(ifnms,3)==1
extrapol = cell(1,4);
extrapol{1} = missedSettingsets{ceil(str2double(name(11:end-14))/3)};
fignr = get(gcf,'number')+1;
figure(fignr);
end
imodel = mod(str2double(name(11:end-14))-1,3)+1;
ydata = squeeze(requiredAmp(1,:,:))*1e9;
xdata = logspace(0,floor(length(ydata)/2),length(ydata))';
[fitresult,coeff,gof] = FitPower(xdata,ydata);
coeff = round(coeff,2,'significant');
fun = @(x) coeff(1).*x.^coeff(2)+coeff(3);
namefun = ['fit: ',num2str(coeff(1)),' x^{',num2str(coeff(2)),'} + ',num2str(coeff(3))];
subplot(2,2,imodel)
plot(xdata,ydata)
set(gca,'xscale','log')
set(findobj(gca,'type','line'),'LineStyle','none')
set(findobj(gca,'type','line'),'Marker','*')
hold on
plot(logspace(0,floor(length(ydata)/2),101),fun(logspace(0,floor(length(ydata)/2),101)),'DisplayName',namefun)
hold off
ylabel('threshold amplitude [nm]')
xlabel('number of dipoles')
title(Titles{imodel});
legend('show')
if imodel == 1
    legend({'RMS threshold = 10%',namefun})
end
if imodel == 3
    extrapol{1+imodel} = fun(10^4);
    set(gcf,'position',[39,3,1904,963]);
    mtit(extrapol{1})
    T = [T;extrapol];
else
    extrapol{1+imodel} = fun(10^6);
end    
end
cd(OLDdirectory);
%%
set(findobj('type','axes'),'fontsize',16)
if isfile('ResultsStrngthcurve.csv')
     overwrite_flag = input('overwrite Results file? 1/0 \n');
     if overwrite_flag
         writetable(T,'ResultsStrngthcurve.csv')
     end
else
     writetable(T,'ResultsStrngthcurve.csv')
end
