function calcStrengthCurve(Input)
% Use this function for HPC measurements
disp('start function')
if ischar(Input)    
Input_temp = load([Input,'.mat']);
Input = Input_temp.Input;
end
Amp_Input = Input.Amp;Totaldps = Input.Totaldps;dps_run = Input.dps_run; Settings = Input.Settings; Fignr = Input.Fignr; SettingsStr = Input.SettingsStr;
HPC_flag = Input.HPC_flag; Models = Input.Models; Solutiontypes = Input.Solutiontypes; RPOI = Input.RPOI;
if HPC_flag
    try
        c = parcluster('local');
        pool = parpool(c.NumWorkers);
    end
end
txtbxstr = gentxtbxstr(Settings);
%Models={'Human scalp','Human air','Human sub scalp','Human cortical','Mouse scalp','Mouse air'};
%Solutiontypes = {'4SphereS8.7R25','4SphereS8.7R25','4SphereS8.7R25','4SphereS8.7R25','Mouse4Sphere~fair0','Mouse4Sphere~fair1'};
%RPOI = [0.082,0.087,0.075,0.07,0.0059,0.0069];
% Models={'Human scalp','Human cortical','Mouse scalp'};
% Solutiontypes = {'4SphereS8.7R25','4SphereS8.7R25','Mouse4Sphere~fair0'};
% RPOI = [0.082,0.07,0.0059];
maxIter = 10;
%Thresholds = [0.1,0.05];
Thresholds = [0.1];
maxosc = Amp_Input(end);
totaliters = length(Thresholds)*length(Totaldps)*maxIter;
try 
    Models{1};
catch
    Models = cellstr(Models);
end
try 
    Solutiontypes{1};
catch
    Solutiontypes = cellstr(Solutiontypes);
end

DEBUG = 0;
PLOT_flag = 0;
LOCALRUN = 0;
%% Debug mode
if DEBUG
    disp('DEBUG modus: limited number of dipoles !!!!')
    disp('')
    Totaldps = Totaldps(7:8);
    dps_run = dps_run(7:8);
    %Amp_Input = Amp_Input([1,end]);
    HPC_flag = 0;
    PLOT_flag = 0;
elseif LOCALRUN
    disp('LOCAL RUN!!!!!')
    disp('___________________________')
    Amp_Input = logspace(0,5,6).*1e-9;
    Totaldps = round(logspace(0,4,5));%round(logspace(0,5+double(isettings>=9 & isettings<=12),6+2*+double(isettings>=9 & isettings<=12)));
    orders = -floor(log10(Totaldps));
    for i=1:length(Totaldps)
        Totaldps(i) = round(Totaldps(i),orders(i));
    end
    dps_run = min(1e3,Totaldps);
end

for imodel = 1:length(Models)
    disp(['Startingmodel: ',Models{imodel}]);
    if contains(lower(Models{imodel}),'mouse')
        Amp = horzcat(1e-10,Amp_Input);
    else
        Amp = Amp_Input; 
    end
    [X,Y] = meshgrid(Totaldps,Amp);
    Y = Y*1e9;
    %fix rng
    rngseed = rng('shuffle');
for iAmp = 1:length(Amp)
    for idps = 1:length(Totaldps)
        [RMS,SNR] = investBiologicalNoise(horzcat({'Totaldps',Totaldps(idps),'dps_run',dps_run(idps),...
            'Aus',Amp(iAmp),'ParallelCompute',dps_run(idps)>100,'SolutionType',Solutiontypes{imodel},...
            'POIs',RPOI(imodel)*[0,0,1],'scale',0,'PLOT',PLOT_flag},Settings));
        RMSnMAT(imodel,iAmp,idps) = RMS(1).RMSnorm; RMSoscMAT(imodel,iAmp,idps) = RMS(1).RMSosc;
        SNRMAT(imodel,iAmp,idps) = SNR(2,1);
    end
end
%% Plot results
figure(Fignr+floor((imodel-1)/3))
subplot(min(3,length(Models)),3,mod(imodel-1,3)*3+1)
surf(X,Y,squeeze(RMSnMAT(imodel,:,:)),'EdgeColor','none','FaceColor','interp')
set(gca,'xscale','log'); set(gca,'yscale','log');
xlabel('# dipoles')
xlim([min(Totaldps),max(Totaldps)])
ylabel({['\bf', Models{imodel}],'','Vibration amplitude [nm]'})
colormap('jet')
colorbar
view(0,90)
title('RMS error')
subplot(min(3,length(Models)),3,2+mod(imodel-1,3)*3)
surf(X,Y,squeeze(RMSoscMAT(imodel,:,:)),'EdgeColor','none','FaceColor','interp')
set(gca,'xscale','log'); set(gca,'yscale','log');
xlabel('# dipoles')
xlim([min(Totaldps),max(Totaldps)])
ylabel('Vibration amplitude [nm]')
colormap('jet')
colorbar
view(0,90)
title('RMS error')
subplot(min(3,length(Models)),3,3+mod(imodel-1,3)*3)
surf(X,Y,squeeze(SNRMAT(imodel,:,:)),'EdgeColor','none','FaceColor','interp')
set(gca,'xscale','log'); set(gca,'yscale','log');
xlabel('# dipoles')
ylabel('Vibration amplitude [nm]')
minlim = Totaldps(Totaldps>1);
xlim([minlim(1),max(Totaldps)])
colormap('jet')
colorbar
view(0,90)
title('SNR')
set(gcf,'position',[-1919,41,1920,963])
if ~mod(imodel,3) || imodel == length(Models)
mtit(SettingsStr,'xoff',-0.01,'yoff',.05,'fontsize',14)
end
%%
    figure(Fignr+floor((length(Models)-1)/3)+1)
    sp{imodel} = subplot(ceil(length(Models)/2),1+double(length(Models)>1),imodel);
    h = waitbar(0,'start');
    hiter = 0;
for ithresh=1:length(Thresholds)
    disp(['start finding threshold: ',num2str(Thresholds(ithresh))]);
    RMSOI = Thresholds(ithresh);
for idps = 1:length(Totaldps)
    superval = Amp(RMSnMAT(imodel,:,idps)<=RMSOI);
    iter = 0;
    findhigherampiter = 0;
    if isempty(superval)
        superval = maxosc;
        newval = superval;
        findhigheramp = 1;        
    else
        superval = superval(1);
        findhigheramp = 0;
    end
    lowerval = Amp(RMSnMAT(imodel,:,idps)>=RMSOI);
    if isempty(lowerval)
        lowerval = superval;
    else
        lowerval = lowerval(end);
    end
    flag = abs(superval-lowerval)>1e-9;
    while flag && iter<maxIter && findhigherampiter<10
        % reset random generator back to original seed => threshold search not
        % dependend on random generation each iteration
        rng(rngseed)
        if findhigheramp
            newval = newval*10
            iter = 0;
            findhigherampiter = findhigherampiter +1;
        else
            newval = 10^((log10(superval)+log10(lowerval))/2)
        end
        [RMS,SNR] = investBiologicalNoise(horzcat({'POIs',RPOI(imodel)*[0,0,1],'Totaldps',Totaldps(idps),...
            'dps_run',dps_run(idps),'Aus',newval,'ParallelCompute',dps_run(idps)>100,...
            'SolutionType',Solutiontypes{imodel},'scale',0},Settings));
        RMSn = RMS(1).RMSnorm;
        if RMSn<=RMSOI
            superval = newval;
            findhigheramp = 0;
        else
            lowerval = newval;
        end
        flag = abs(superval-lowerval)>1e-9;
        waitbar((hiter*maxIter+iter)/totaliters,h, [num2str((hiter*maxIter+iter)/totaliters*100),'%'])
        iter = iter+1;
    end
    hiter = hiter+1;
    requiredAmp(imodel,ithresh,idps) = superval;
    requiredAmplower(imodel,ithresh,idps) = lowerval;
    itermat(imodel,ithresh,idps) = iter;
    %save intermediate results
    idxdoublepoint = find(SettingsStr==':');
    if isempty(idxdoublepoint)
        save_str_int = SettingsStr;
    else
        save_str_int = SettingsStr([1:idxdoublepoint-1,idxdoublepoint+2:end]);
    end
    save(['Results_',save_str_int,'_',datestr(now,'mm-dd-yy_HHMM'),'.mat'],'requiredAmp','requiredAmplower','itermat');
end
plot(Totaldps,squeeze(requiredAmp(imodel,ithresh,:)*1e9))
hold on
end
hold off
set(gca,'xscale','log');
xlabel('number of dipoles')
ylabel('threshold amplitude [nm]')
title(Models{imodel})
if mod(imodel+1,2) || imodel == length(Models)
legend(arrayfun(@(x) ['RMS threshold = ',num2str(x*100),'%'],Thresholds,'UniformOutput',false))
end
close(h)
end
set(gcf,'position',[-1919,41,1920,963])
pause(0.1)
pos = get(sp{1}, 'position');
dim = pos.*[1 1+0.5*pos(4)/pos(2) 0.5 0.5];
annotation(gcf, 'textbox', dim, 'String', txtbxstr,'FitBoxToText','on')
mtit(SettingsStr,'xoff',-0.001,'yoff',.05,'fontsize',14)
disp('');
disp(['ended, start saving']);
if HPC_flag
    idxdoublepoint = find(SettingsStr==':');
    if isempty(idxdoublepoint)
        save_str_int = SettingsStr;
    else
        save_str_int = SettingsStr([1:idxdoublepoint-1,idxdoublepoint+2:end]);
    end
    save(['Results_',save_str_int,'_',datestr(now,'mm-dd-yy_HHMM'),'.mat'],'requiredAmp','requiredAmplower','itermat');
    if floor((length(Models)-1)/3)+1==2
        Figtitles = {'Surfplots','Surfplots2','Strngthcurves'};
    else
        Figtitles = {'Surfplots','Strngthcurves'};
    end
    for ifigsave=0:floor((length(Models)-1)/3)+1
        figstr = [Figtitles{mod(ifigsave,3)+1},save_str_int,'_',datestr(now,'mm-dd-yy_HHMM'),'.fig'];
        figure(Fignr+ifigsave)
        pause(0.01)
        savefig(gcf,figstr)
        pause(0.01)
    end
end
    close all
end
