% do pulsed study
clear all

Input = load('./Inputs/threshold05/Input_cEG_Amp_fkApdpset_NCfixed_dps16_nos4l_v2_1.mat');
Input = Input.Input;
Input.kAus = Input.kAus(1:4);
Input.Totaldps = Input.Totaldps(1:4);
Input.dps_run = Input.dps_run(1:4);
Input.Settings{14} = [Input.Settings{14},'3']; 
pulsed.flag = 1;
prfs = [2000,4000,8000,10000];
dcs = [0.01,0.04,0.05,0.08,0.1,0.12,0.15,0.16,0.2,1];
Input.Amp = Input.Amp(1:2:end);
Input.save_flag = 0;

Results = {};
for iprf =1:length(prfs)
    for idc = 1:length(dcs)
        pulsed.prp = 1/prfs(iprf);
        pulsed.dc = dcs(idc);
        if iprf==1 && idc == 1
            Input.Settings = horzcat(Input.Settings,{'pulsed',pulsed});
        else
            Input.Settings{end-1} = pulsed;
        end
        [RMSall,SNRall,Outall,TSTOPall,Totaldps,requiredAmp,requiredAmplower,itermat,cSC_mode] = calcErrorGrid(Input);
        Results{end+1} = {RMSall,SNRall,Outall,TSTOPall,Totaldps,requiredAmp,requiredAmplower,itermat,cSC_mode};
    end
end
%%
reQAmp = [];
for iprf =1:length(prfs)
    for idc = 1:length(dcs)
        reQAmp(iprf,idc,:) = Results{(iprf-1)*length(dcs)+idc}{6};
    end
end
figure()
for iprf =1:5
    subplot(2,3,iprf)
    imagesc(log10(squeeze(reQAmp(:,:,iprf))))
    colorbar
    %caxis(log10([min(reQAmp(:)),max(reQAmp(:))]))
    title(num2str(num2str(iprf)))
end

%%
close all
myfun =@(f,Tp,Ts) Tp*sin(pi*f*Tp)./(Ts*f*Tp);
fs = logspace(5,7,11);
figure
for iprf =1:length(prfs)
subplot(2,2,iprf)
for idc = 1:length(dcs)
    Tp = dcs(idc)/prfs(iprf)
    Ts = 1/prfs(iprf)
plot(fs,myfun(fs,Tp,Ts))
hold on
set(gca,'xscale','log')
end

end
