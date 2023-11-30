oldFolder = cd('Figures biological noise');
fignames = {'Currents','Sphere','OverviewMeasuredSignal','','','','ReconstructedSignals','SpectrumPOI1','SpectrumPOI2','SpectrumPOI3'};
nrdipoles = '100000RadialDipoles100nm_';
mkdir(nrdipoles(1:end-1))
cd(nrdipoles(1:end-1))
for i = [1,2,3,7,8,9,10]
    figure(i)
    set(gcf,'position',[-1919,41,1920,963])
    savefig([nrdipoles,fignames{i}])
end
cd(oldFolder)