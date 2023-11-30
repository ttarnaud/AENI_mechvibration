clear all
close all
clc

surfplot_flag = 1;
bargraph_flag = 0;
lineplot_flag = 1;
load('./Results/SCs_all_cEG_Amp_nos4l_fkApdpset_v1.mat')

Totaldps = Inputall(1).Input.Totaldps;
Amp = Inputall(1).Input.Amp;
Amp = Amp*1e6; %amp in um

types = {'-','--',':','-.'};
Markers = {'o','x','d','s'};

iM = size(requiredAmpAll,1);
Models = {};
ROIs = {};
nLayers = [];
for i=1:iM
    Models = horzcat(Models,Inputall(i).Input.Model);
    ROIs = horzcat(ROIs,Inputall(i).Input.ROI);
    idx = find(strcmpi(Inputall(i).Input.Settings,'nLayers'));
    nLayers = horzcat(nLayers,Inputall(i).Input.Settings{idx+1});
end

uModels = unique(Models);
uROI = unique(ROIs);
unLayers = unique(nLayers);
Colors = lines(length(uModels));

CM = viridis; %CM = flip(CM);
RMS_str = {'RMSDOIvrpOSCvr'}
models_idx = 1:length(RMSALL);
if surfplot_flag
    for i=models_idx
        for iRMS = 1:length(RMS_str)
            X = [];
            Y = [];
            Z2 = [];
            Zmax = [];
            ZmeanPSO = [];
            Zmean = [];
            for iAmp = 1:length(Amp)
                for idps = 1:length(Totaldps)
                    Y(iAmp,idps) = Amp(iAmp);
                    X(iAmp,idps) = Totaldps(idps);
                    Z2(iAmp,idps) = RMSALL(i).RMS(iAmp,idps).RMS(2).(RMS_str{iRMS});
                    Zmax(iAmp,idps) = max([RMSALL(i).RMS(iAmp,idps).RMS(:).(RMS_str{iRMS})]);
                    ZmeanPSO(iAmp,idps) = max([RMSALL(i).RMS(iAmp,idps).RMS(7).(RMS_str{iRMS})]);
                    Zmean(iAmp,idps) = max([RMSALL(i).RMS(iAmp,idps).RMS(6).(RMS_str{iRMS})]);
                end
            end
            figure
            
            colormap(CM)
            surf_set = {'facecolor','interp'};
            sp_set = {{'xscale','yscale','xlim','ylim'},{'log','log',[min(X(:)),max(X(:))],[min(Y(:)),max(Y(:))]}};
            nlayers = Inputall(i).Input.Settings{find(strcmpi(Inputall(i).Input.Settings,'nLayers'))+1};
            master_title = sprintf('%s - %s, nlayers = %i: %s ',Inputall(i).Input.Model,Inputall(i).Input.ROI,nlayers,RMS_str{iRMS});
            
            subplot(2,2,1)
            surf(X,Y,Z2,Z2,surf_set{:})
            hold on
            [x,y,z] = limline(X,Y,Z2,0.1);
            plot3(x,y,z+1,'r')
            hold off
            view(0,90)
            set(gca,sp_set{1},sp_set{2})
            colorbar
            xlabel('# dipoles')
            ylabel('Vibration amplitude [\mum]')
            title('RMSE at POI on z-axis')
            
            subplot(2,2,2)
            surf(X,Y,Zmax,Zmax,surf_set{:})
            hold on
            [x,y,z] = limline(X,Y,Zmax,0.1);
            plot3(x,y,z+1,'r')
            hold off
            view(0,90)
            set(gca,sp_set{1},sp_set{2})
            colorbar
            xlabel('# dipoles')
            ylabel('Vibration amplitude [\mum]')
            title('max RMSE')
            
            subplot(2,2,3)
            surf(X,Y,Zmean,Zmean,surf_set{:})
            hold on
            [x,y,z] = limline(X,Y,Zmean,0.1);
            plot3(x,y,z+1,'r')
            hold off
            view(0,90)
            set(gca,sp_set{1},sp_set{2})
            colorbar
            xlabel('# dipoles')
            ylabel('Vibration amplitude [\mum]')
            title('RMSE mean all POIs ')
            
            subplot(2,2,4)
            surf(X,Y,ZmeanPSO,ZmeanPSO,surf_set{:})
            hold on
            [x,y,z] = limline(X,Y,ZmeanPSO,0.1);
            plot3(x,y,z+1,'r')
            hold off
            view(0,90)
            set(gca,sp_set{1},sp_set{2})
            colorbar
            xlabel('# dipoles')
            ylabel('Vibration amplitude [\mum]')
            title('RMSE mean POIs same O')
            mtit(master_title);
        end
        
    end
end
RMS_str = {'RMS','RMSDOIvrpOSCvr','RMSDOIvrpstatnoise'}
models_idx = 1:length(RMSALL);

if bargraph_flag || lineplot_flag
    Method = 'lin';
    for iRMS = 1:length(RMS_str)
        Val2.(RMS_str{iRMS}) = [];
        ValmPSO.(RMS_str{iRMS}) = [];
        for i=models_idx
            idx_c = find(strcmpi(Models{i},uModels));
            idx_lt = find(strcmpi(ROIs{i},uROI));
            idx_m = find(nLayers(i)==unLayers);
            
            
            X = [];
            Y = [];
            Z2 = [];
            Zmax = [];
            ZmeanPSO = [];
            Zmean = [];
            for iAmp = 1:length(Amp)
                for idps = 1:length(Totaldps)
                    
                    Y(iAmp,idps) = Amp(iAmp);
                    X(iAmp,idps) = Totaldps(idps);
                    Z2(iAmp,idps) = RMSALL(i).RMS(iAmp,idps).RMS(2).(RMS_str{iRMS});
                    Zmax(iAmp,idps) = max([RMSALL(i).RMS(iAmp,idps).RMS(:).(RMS_str{iRMS})]);
                    ZmeanPSO(iAmp,idps) = max([RMSALL(i).RMS(iAmp,idps).RMS(7).(RMS_str{iRMS})]);
                    Zmean(iAmp,idps) = max([RMSALL(i).RMS(iAmp,idps).RMS(6).(RMS_str{iRMS})]);
                end
                
            end
            idx = length(Totaldps);
            zval = 10^interp1(Z2(:,idx),log10(Y(:,idx)),0.1,Method,'extrap');
            Val2.(RMS_str{iRMS})(idx_m,idx_c,idx_lt) = zval;
            zmPSO = 10^interp1(Z2(:,idx),log10(Y(:,idx)),0.1,Method,'extrap');
            ValmPSO.(RMS_str{iRMS})(idx_m,idx_c,idx_lt) = zmPSO;
        end
        if bargraph_flag 
            figure
            subplot(2,2,1)
            bar(Val2.(RMS_str{iRMS})(:,:,1));
            set(gca,{'yscale','xticklabel'},{'log',unLayers})
            title(uROI{1})
            subplot(2,2,2)
            bar(Val2.(RMS_str{iRMS})(:,:,2));
            set(gca,{'yscale','xticklabel'},{'log',unLayers})
            legend(uModels)
            title(uROI{2})
            
            subplot(2,2,3)
            bar(ValmPSO.(RMS_str{iRMS})(:,:,1));
            set(gca,{'yscale','xticklabel'},{'log',unLayers})
            title(uROI{1})
            subplot(2,2,4)
            bar(ValmPSO.(RMS_str{iRMS})(:,:,2));
            set(gca,{'yscale','xticklabel'},{'log',unLayers})
            legend(uModels)
            title(uROI{2})
        end
        if lineplot_flag
            fignr = get(gcf,'number')+1;
            s1_flag = 1;
            s2_flag = 1;
            for i=models_idx
                idx_c = find(strcmpi(Models{i},uModels));
                idx_lt = find(strcmpi(ROIs{i},uROI));
                idx_m = find(nLayers(i)==unLayers);
                
                
                X = [];
                Y = [];
                Z2 = [];
                Zmax = [];
                ZmeanPSO = [];
                Zmean = [];
                for iAmp = 1:length(Amp)
                    for idps = 1:length(Totaldps)
                        
                        Y(iAmp,idps) = Amp(iAmp);
                        X(iAmp,idps) = Totaldps(idps);
                        Z2(iAmp,idps) = RMSALL(i).RMS(iAmp,idps).RMS(2).(RMS_str{iRMS});
                        Zmax(iAmp,idps) = max([RMSALL(i).RMS(iAmp,idps).RMS(:).(RMS_str{iRMS})]);
                        ZmeanPSO(iAmp,idps) = max([RMSALL(i).RMS(iAmp,idps).RMS(7).(RMS_str{iRMS})]);
                        Zmean(iAmp,idps) = max([RMSALL(i).RMS(iAmp,idps).RMS(6).(RMS_str{iRMS})]);
                    end
                    
                end
                for idps=1:length(Totaldps)
                    zval(idps) = 10^interp1(Z2(:,idps),log10(Y(:,idps)),0.1,Method,'extrap');                    
                    zmPSO(idps) = 10^interp1(Z2(:,idps),log10(Y(:,idps)),0.1,Method,'extrap');
                    
  
                    
                end
                ltype = [Markers{idx_m},'-'];
                
                figure(fignr)
                subplot(1,2,idx_lt)
                p = plot(Totaldps,zval,ltype,'color',Colors(idx_c,:),'markersize',10);
                if s1_flag && idx_lt==1
                    hold on;
                    title(ROIs{i})
                    xlabel('nr. dipoles')
                    ylabel('Amp [um]')
                    s1_flag = 0;
                    ylim([10^-3,10^3])
                end
                if s2_flag && idx_lt==2
                    hold on;
                    title(ROIs{i})
                    xlabel('nr. dipoles')
                    s2_flag = 0;
                    ylim([10^-3,10^3])
                end
                set(gca,{'xscale','yscale'},{'log','log'})
                set(get(get(p,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
            end
            for i=1:length(uModels)
                plot(nan,nan,'-','color',Colors(i,:),'DisplayName',[uModels{i}(1:5),'_{',uModels{i}(7:end),'}']);
            end
            for i=1:length(unLayers)
                plot(nan,nan,Markers{i},'color','k','DisplayName',sprintf('# Layers = %i',unLayers(i)));
            end
            hold off
            legend('show','location','best')
            set(findall(gcf,'type','axes'),'fontsize',18)
            
        end
        
    end
end


function [x,y,z] = limline(X,Y,Z,Lim)
idxZ = Z<=Lim;
idxline = idxZ & circshift(~idxZ,-1,2);
x=X(idxline); x=x(:);
y=Y(idxline); y = y(:);
z = Z(idxline); x=x(:);

end