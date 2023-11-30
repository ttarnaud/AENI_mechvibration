
%Collect all info
intm1 = load('./Results/SCs_all_cEG_Amp_nos4l_fkApdpset_v2.mat');
intm2 = load('./Results/SCs_all_cEG_Amp_NCfixed2_nos4l_fkApdpset_v2.mat');
intm3 = load('./Results/SCs_all_cEG_Amp_NCfixed2_sOSC_nos4l_fkApdpset_v2.mat');
intm4 = load('./Results/SCs_all_cEG_Amp_NCfixed3_sOSC_nos4l_fkApdpset_v2.mat');

dps1 = intm1.Inputall(1).Input.Totaldps;
dps2 = intm2.Inputall(1).Input.Totaldps;
dps3 = intm3.Inputall(1).Input.Totaldps;
dps4 = intm4.Inputall(1).Input.Totaldps;

dps = unique([dps1,dps2,dps3,dps4]);

RequiredAmp1 = intm1.requiredAmpAll;
RequiredAmp2 = intm2.requiredAmpAll;
RequiredAmp3 = intm3.requiredAmpAll;
RequiredAmp4 = intm4.requiredAmpAll;

M = max([size(RequiredAmp1,1),size(RequiredAmp1,1),size(RequiredAmp1,1),size(RequiredAmp1,1)]);
RequiredAmp = nan(M,length(dps),4);


nLayers = [3,5,10];
Models = {'human_scalp','human_cortical','human_air','mouse_fair0'};
ROIs = {'Deep','Cortex'};
[nL_mat,ROI_mat,Models_mat] = ndgrid(nLayers,ROIs,Models);
nL_mat = nL_mat(:)'; ROI_mat = ROI_mat(:)'; Models_mat = Models_mat(:)';

[~,idx_dps,idx_dps2] = intersect(dps,dps1);
for i=1:length(intm1.Inputall)
    Settings = intm1.Inputall(i).Input.Settings;
    info1(i).nLayers = Settings{find(strcmpi(Settings,'nLayers'))+1};
    info1(i).Model = intm1.Inputall(i).Input.Model;
    info1(i).ROI = intm1.Inputall(i).Input.ROI;
    
    
    idx = info1(i).nLayers==nL_mat & strcmpi(info1(i).Model,Models_mat) &...
        strcmpi(info1(i).ROI,ROI_mat);
    RequiredAmp(idx,idx_dps,1) = RequiredAmp1(i,idx_dps2);
end

[~,idx_dps,idx_dps2] = intersect(dps,dps2);
for i=1:length(intm2.Inputall)
    Settings = intm2.Inputall(i).Input.Settings;
    info2(i).nLayers = Settings{find(strcmpi(Settings,'nLayers'))+1};
    info2(i).Model = intm2.Inputall(i).Input.Model;
    info2(i).ROI = intm2.Inputall(i).Input.ROI;
    
    idx = info2(i).nLayers==nL_mat & strcmpi(info2(i).Model,Models_mat) & strcmpi(info2(i).ROI,ROI_mat);
    RequiredAmp(idx,idx_dps,2) = RequiredAmp2(i,idx_dps2);
end

[~,idx_dps,idx_dps2] = intersect(dps,dps3);
for i=1:length(intm3.Inputall)
    Settings = intm3.Inputall(i).Input.Settings;
    info3(i).nLayers = Settings{find(strcmpi(Settings,'nLayers'))+1};
    info3(i).Model = intm3.Inputall(i).Input.Model;
    info3(i).ROI = intm3.Inputall(i).Input.ROI;
    
    idx = info3(i).nLayers==nL_mat & strcmpi(info3(i).Model,Models_mat) & strcmp(info3(i).ROI,ROI_mat);
    RequiredAmp(idx,idx_dps,3) = RequiredAmp3(i,idx_dps2);
    
end

[~,idx_dps,idx_dps2] = intersect(dps,dps4);
for i=1:length(intm4.Inputall)
    Settings = intm4.Inputall(i).Input.Settings;
    info4(i).nLayers = Settings{find(strcmpi(Settings,'nLayers'))+1};
    info4(i).Model = intm4.Inputall(i).Input.Model;
    info4(i).ROI = intm4.Inputall(i).Input.ROI;
    
    idx = info4(i).nLayers==nL_mat & strcmpi(info4(i).Model,Models_mat) & strcmpi(info4(i).ROI,ROI_mat);
    RequiredAmp(idx,idx_dps,4) = RequiredAmp4(i,idx_dps2);
end

%%
row1 = {'\alpha_2 = 2','3','3','4'};
row2 = {'PSD(1MHZ) = 1e-11','1e-14','1e-14','1e-17'};
xvals = [1,3,4]
labelArray = [row1(xvals);row2(xvals)];
labelArray = strjust(pad(labelArray),'center');
tickLabels = strtrim(sprintf('%s\\newline%s\n', labelArray{:}))

types = {'-','--',':','-.'};
Markers = {'o','x','d','s'};
Colors = lines(length(Models));

dps_nr = [3,5,7,9];
YVAL = [];
Rsqval = [];
alphaval_y50nm = [];
Model_rval = [];
PSDval_y50nm = [];
figure
for idps = 1:length(dps_nr)
    for iM = 1:length(Models_mat)
        
        
    idx_c = find(strcmpi(Models_mat{iM},Models));
    idx_lt = find(strcmpi(ROI_mat{iM},ROIs));
    idx_m = find(nL_mat(iM)==nLayers);
    ltype = [Markers{idx_m},types{idx_lt}];
    yval = squeeze(RequiredAmp(iM,dps_nr(idps),xvals));
    yval = yval*1e6;
    if length(yval(~isnan(yval)))>1
        yintm = log10(yval'/mean(yval));
        YVAL = [YVAL;yintm];
        Model_rval = [Model_rval;[idps,iM]];
        b = [ones(numel(yval),1),[1;2;3]]\yintm';
        Rval = 1-sum((yintm'-[ones(numel(yval),1),[1;2;3]]*b).^2)/sum((yintm'-mean(yintm)).^2);
        Rsqval = [Rsqval;Rval];
        % find x for y =100nm
        yintm = log10(yval');
        X0 = [ones(numel(yval),1),[2;3;4]];
        b = X0\yintm';
        xval_y50nm = (log10(0.050)-b(1))/b(2);
        PSDval = calcPSD1MHz(xval_y50nm);
        alphaval_y50nm = [alphaval_y50nm;xval_y50nm];
        PSDval_y50nm = [PSDval_y50nm;PSDval];
    end
    subplot(2,2,idps)
    p = plot(1:length(xvals),log10(yval),ltype,'color',Colors(idx_c,:),'markersize',10);
    if iM==1; hold on; end
    
    %set(gca,{'xscale','yscale'},{'log','log'})
    %set(gca,{'yscale'},{'log'});
    set(get(get(p,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
    ax = gca;
    ax.XLim = [0,length(xvals)+1];
    ax.XTick = 1:length(xvals);
    ax.XTickLabel = tickLabels;
    
    end
    title(['nr. dipoles =',  num2str(dps(dps_nr(idps)))])
end
subplot(2,2,4)
for i=1:length(Models)
    plot(nan,nan,'-','color',Colors(i,:),'DisplayName',[Models{i}(1:5),'_{',Models{i}(7:end),'}']);
end
for i=1:length(ROIs)
    plot(nan,nan,types{i},'color','k','DisplayName',ROIs{i});
end
for i=1:length(nLayers)
    plot(nan,nan,Markers{i},'color','k','DisplayName',sprintf('# Layers = %i',nLayers(i)));
end
hold off
legend('show')
subplot(2,2,1)
ylabel('Vibration Amplitude (log) [\mum]')
subplot(2,2,3)
ylabel('Vibration Amplitude (log) [\mum]')

% linear regresion
XVAL = repmat([1,2,3],size(YVAL,1),1);
b = [ones(numel(XVAL),1),XVAL(:)]\YVAL(:)

Rsq = 1-sum((YVAL(:)-[ones(numel(XVAL),1),XVAL(:)]*b).^2)/sum((YVAL(:)-mean(YVAL(:))).^2)
Rsq_min = min(Rsqval)
Rsq_max = max(Rsqval)
figure
scatter(XVAL(:),YVAL(:));
hold on
plot([1,2,3],[ones(3,1),[1,2,3]']*b,'k-')
hold off
xlim([0,4])


function PSD = calcPSD1MHz(alpha2)
%noise not included
wn = 14;
theta = 0.6;         %max is reached between wn and theta*wn
w2 = 10^3;                      %breakpoint higher powerlaw
alpha1 = 2;
%noise2 =@(x) 10.^(randn(1,length(x)));  %different noise methods
%were tested see script An_spec
noise3 = @(x) lognrnd(0,0.5,1,length(x));
y = @(x) (1+x/(theta*wn))./((1+(x/wn).^alpha1).*(1+(x/w2).^alpha2));
fun = @(x) y(x).*noise3(x);
PSD = y(1e6);
end