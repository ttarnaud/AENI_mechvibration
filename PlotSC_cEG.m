%% load Strenght cruves and plot of simulations Input_cEG_nos4l_v1_1-24
load('./Results/SCs_all_cEG_Amp_nos4l_fkApdpset_v1.mat')

Totaldps = Inputall(1).Input.Totaldps;
Amp = Inputall(1).Input.Amp;

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



figure
for i=1:iM
    idx_c = find(strcmpi(Models{i},uModels));
    idx_lt = find(strcmpi(ROIs{i},uROI));
    idx_m = find(nLayers(i)==unLayers);
    ltype = [Markers{idx_m},types{idx_lt}];
    yval = requiredAmpAll(i,:);
    yval = yval*1e6;
    
    [xnew,ynew] = exclude_vals(Totaldps,yval,1e4,1e3);
    
    p = plot(xnew,ynew,ltype,'color',Colors(idx_c,:),'markersize',10)
    if i==1; hold on; end
    
    set(gca,{'xscale','yscale'},{'log','log'})
    set(get(get(p,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
    
    
end

for i=1:length(uModels)
    plot(nan,nan,'-','color',Colors(i,:),'DisplayName',[uModels{i}(1:5),'_{',uModels{i}(7:end),'}']);
end
for i=1:length(uROI)
    plot(nan,nan,types{i},'color','k','DisplayName',uROI{i});
end
for i=1:length(unLayers)
    plot(nan,nan,Markers{i},'color','k','DisplayName',sprintf('# Layers = %i',unLayers(i)));
end
hold off
legend('show','location','best')

s1_flag = 1;
s2_flag = 1;
figure
for i=1:iM
    idx_c = find(strcmpi(Models{i},uModels));
    idx_lt = find(strcmpi(ROIs{i},uROI));
    idx_m = find(nLayers(i)==unLayers);
    ltype = [Markers{idx_m},'-'];
    yval = requiredAmpAll(i,:);
    yval = yval*1e6;
    subplot(1,2,idx_lt)
    [xnew,ynew] = exclude_vals(Totaldps,yval,1e4,1e3);
    p = plot(xnew,ynew,ltype,'color',Colors(idx_c,:),'markersize',10)
    if s1_flag && idx_lt==1
    hold on;
    title(ROIs{i})
    xlabel('nr. dipoles')
    ylabel('vibration Amp [um]')
    end
    if s2_flag && idx_lt==2
    hold on;
    title(ROIs{i})
    xlabel('nr. dipoles')
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


% 8 subplots
fignr = get(gcf,'number');
figure
for i=1:iM
    idx_c = find(strcmpi(Models{i},uModels));
    idx_lt = find(strcmpi(ROIs{i},uROI));
    idx_m = find(nLayers(i)==unLayers);
    figure(fignr+idx_m)
    ltype = [Markers{1},'-'];    
    yval = requiredAmpAll(i,:);
    yval = yval*1e6;
    subplot(1,2,idx_lt)
    [xnew,ynew] = exclude_vals(Totaldps,yval,1e4,1e3);
    p = plot(xnew,ynew,ltype,'color',Colors(idx_c,:),'markersize',10)
    if s1_flag && idx_lt==1
    hold on;
    title(ROIs{i})
    xlabel('nr. dipoles')
    ylabel('vibration Amp [um]')
    end
    if s2_flag && idx_lt==2
    hold on;
    title(ROIs{i})
    xlabel('nr. dipoles')
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



function [xnew,ynew] = exclude_vals(x,y,limitx,limity)

idx_x = x<=limitx;
x_intm = x(idx_x);
y_intm = y(idx_x);

idx_y = y_intm<=limity;
ynew = [y_intm(idx_y),y(~idx_x)];
xnew = [x_intm(idx_y),x(~idx_x)];


end