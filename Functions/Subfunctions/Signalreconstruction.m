function [RMS,SNR,timeReconS_f,reconsSignalVR_f,reconsSignalVRDOI_f,reconsSignalVRdv_f,...
    reconsSignalVRpthn_f,reconsSignalVRdvpov_f,reconsSignalVRdvpsn_f,VR_f,VRpthn_f,VRDOI_f,VRDOIvibr_f] =...
    Signalreconstruction(VR,VRDOI,VRds,VROSC,VRos,VRsn,ThermalNoise,theta_global,Tsim,POIs,POIsidx,fus,fbandwidth,inputSignal,mirror,PtP,PLOT,HPC_flag,Display,PLOT_RMS)
% in this function the a reconstruction is made of the signal near the fus
% +-fbandwidth. RMS wrt input signal for reconstruction of
% VR,VRDOI,VRDOIvibr (latter should be perfect,second small error due to
% static noise,VR larger errors due to all noise sources) also VRptn VRplus
% thermalnoise, VRdvpthn
%SNR: VRDOIvibr/(VROSCvibr+VRstatnoise)   Signal power vs all noise
%   : VRDOIvibr/(VROSCvibr)               Signal power vs vibr noise
%   :VRDOIvibr/(VRstatnoise)  singal vs tat noise
%   : VRDOIvibr/(VROSCvibr+VRstatnoise+VRptherm)
%VR = VRstatnoise+VROSCvibr+VRDOIvibr
%VRstatnoise = VRstatstat+VROSCstat+VRDOIstat
%VROSC = VROSCvibr+VROSCstat
%VRDOI = VRDOIvibr+VRDOIstat
%VRstatnoise = VRsn, VROSCvibr = VRov, VRDOIvibr = VRdv;
%VRstatstat = VRss; VROSCstat = VRos; VRDOIstat=VRds

%OUTPUT: RMS:structure conatining results of RMS
%        SNr: structure containing results of SNR
%        TimeReconS_f: cell array with time points of reconstructed signal
%        also additional meanPOIs and mean PsO: _f applies to all other
%        outpurs as indication signal of additionals are added

%individual components
VRov = VROSC-VRos;         %Vibrational component of oscillators except DOI
VRdv = VRDOI-VRds;         %vibrational component of DOI
VRpthn = VR+ThermalNoise;  %signal plus Thermalnoise
VRan = VRsn+VRov;          %all noise sources 
VRanpthn = VRan+ThermalNoise; %all noise plus thermal noise
VRdvpthn = VRdv+ThermalNoise;

%Interesting RMS metrics could be RMS DOIvr+OSCvr (result of only
%vibrational effects, case signal static components are negligable)
%DOIvr+Statnoise: case of perfect field
VRdvpov = VRdv+VRov;
VRdvpsn = VRdv+VRsn;




%reconstruct based on complete signal    
[reconsSignalVR,~,timeReconS] = calcFourierGetsignal(VR,Tsim,POIsidx,fus,fbandwidth,inputSignal,mirror,theta_global,PLOT,'POIstoPLOT',PtP);
if ~HPC_flag && PLOT
    mtit('based on all dipoles')
    pause(0.1)
end
%reconstruct based on complete signal plus thermal   
[reconsSignalVRpthn] = calcFourierGetsignal(VRpthn,Tsim,POIsidx,fus,fbandwidth,inputSignal,mirror,theta_global,PLOT,'POIstoPLOT',PtP);
if ~HPC_flag && PLOT
    mtit('based on signal all dipoles + thermal noise')
    pause(0.1)
end
% reconstruct based on only DOI
[reconsSignalVRDOI] = calcFourierGetsignal(VRDOI,Tsim,POIsidx,fus,fbandwidth,inputSignal,mirror,theta_global,0);
% reconstruct based on only DOIvibr
[reconsSignalVRdv,powerSignalVRdv] = calcFourierGetsignal(VRdv,Tsim,POIsidx,fus,fbandwidth,inputSignal,mirror,theta_global,0);
% reconstruct based on only DOIvibrptherm
[reconsSignalVRdvpthn] = calcFourierGetsignal(VRdvpthn,Tsim,POIsidx,fus,fbandwidth,inputSignal,mirror,theta_global,0);
% reconstruct based on all noise
[reconsSignalVRan,powerSignalVRan] = calcFourierGetsignal(VRan,Tsim,POIsidx,fus,fbandwidth,inputSignal,mirror,theta_global,0);
% reconstruct based on all noise pthn
[reconsSignalVRanpthn,powerSignalVRanpthn] = calcFourierGetsignal(VRanpthn,Tsim,POIsidx,fus,fbandwidth,inputSignal,mirror,theta_global,0);
% reconstruct based on static noise
[reconsSignalVRsn,powerSignalVRsn] = calcFourierGetsignal(VRsn,Tsim,POIsidx,fus,fbandwidth,inputSignal,mirror,theta_global,0);
% reconstruct based on OSCvibr
[reconsSignalVRov,powerSignalVRov] = calcFourierGetsignal(VRov,Tsim,POIsidx,fus,fbandwidth,inputSignal,mirror,theta_global,0);
% reconstruct based on DOivr + OSCvr
[reconsSignalVRdvpov] = calcFourierGetsignal(VRdvpov,Tsim,POIsidx,fus,fbandwidth,inputSignal,mirror,theta_global,0);
% reconstruct based on DOivr + statnoise
[reconsSignalVRdvpsn] = calcFourierGetsignal(VRdvpsn,Tsim,POIsidx,fus,fbandwidth,inputSignal,mirror,theta_global,0);


RMS = struct;
SNR = struct;
PowerDOIvr = zeros(1,size(POIs,1)+2);       %Power vibrational component DOI
PowerNoiseAll = zeros(1,size(POIs,1)+2);    % power all noise components
PowerNoiseAllpthn = zeros(1,size(POIs,1)+2); %power all noise plus thermal noise
PowerNoiseOSCvr = zeros(1,size(POIs,1)+2);     %power of vibrational interference
PowerNoiseStatic = zeros(1,size(POIs,1)+2);   %power of all static noise

NPOIs = length(reconsSignalVR);
for iPOI=1:NPOIs
    % SNR
    PowerDOIvr(iPOI) = calcPower(reconsSignalVRdv{iPOI});
    PowerNoiseAll(iPOI) = calcPower(reconsSignalVRan{iPOI});
    PowerNoiseAllpthn(iPOI) = calcPower(reconsSignalVRanpthn{iPOI});
    PowerNoiseOSCvr(iPOI) = calcPower(reconsSignalVRov{iPOI});
    PowerNoiseStatic(iPOI) = calcPower(reconsSignalVRsn{iPOI});
   
    SNR.DOIvr_NoiseAll(1,iPOI) = 10*log10(PowerDOIvr(iPOI)/PowerNoiseAll(iPOI));
    SNR.DOIvr_NoiseAll(2,iPOI) = 10*log10(powerSignalVRdv(iPOI)/powerSignalVRan(iPOI));
    SNR.DOIvr_NoiseAllptherm(1,iPOI) = 10*log10(PowerDOIvr(iPOI)/PowerNoiseAllpthn(iPOI));
    SNR.DOIvr_NoiseAllptherm(2,iPOI) = 10*log10(powerSignalVRdv(iPOI)/powerSignalVRanpthn(iPOI));
    SNR.DOIvr_NoiseOSCvr(1,iPOI) = 10*log10(PowerDOIvr(iPOI)/PowerNoiseOSCvr(iPOI));
    SNR.DOIvr_NoiseOSCvr(2,iPOI) = 10*log10(powerSignalVRdv(iPOI)/powerSignalVRov(iPOI));
    SNR.DOIvr_NoiseStatic(1,iPOI) = 10*log10(PowerDOIvr(iPOI)/PowerNoiseStatic(iPOI));
    SNR.DOIvr_NoiseStatic(2,iPOI) = 10*log10(powerSignalVRdv(iPOI)/powerSignalVRsn(iPOI));
    
    %RMS
    [RMS(iPOI).RMS,~,RMS(iPOI).Q2] = ...
        calcRMS(reconsSignalVR{iPOI},inputSignal,timeReconS{iPOI},mirror(iPOI),PLOT_RMS);
    if ~HPC_flag && PLOT_RMS
        mtit(sprintf('based on all dipoles POI %i',iPOI))
        pause(0.1)
    end
    [RMS(iPOI).RMSptherm,~,RMS(iPOI).Q2ptherm] = ...
        calcRMS(reconsSignalVRpthn{iPOI},inputSignal,timeReconS{iPOI},mirror(iPOI),PLOT_RMS);
    if ~HPC_flag && PLOT_RMS
        mtit(sprintf('based on all dipoles + thermal noise POI %i',iPOI))
        pause(0.1)
    end
    
    [RMS(iPOI).RMSDOI,~,RMS(iPOI).Q2DOI] = ...
        calcRMS(reconsSignalVRDOI{iPOI},inputSignal,timeReconS{iPOI},mirror(iPOI),PLOT_RMS);
    if ~HPC_flag && PLOT_RMS
        mtit(sprintf('based on DOI POI %i',iPOI))
        pause(0.1)
    end
    [RMS(iPOI).RMSDOIvr,~,RMS(iPOI).Q2DOIvr] = ...
        calcRMS(reconsSignalVRdv{iPOI},inputSignal,timeReconS{iPOI},mirror(iPOI),PLOT_RMS);
    if ~HPC_flag && PLOT_RMS
        mtit(sprintf('based vibrational component DOI  POI %i',iPOI))
        pause(0.1)
    end
    [RMS(iPOI).RMSDOIvrpthn,~,RMS(iPOI).Q2DOIvrpthn] = ...
        calcRMS(reconsSignalVRdvpthn{iPOI},inputSignal,timeReconS{iPOI},mirror(iPOI),PLOT_RMS);
    if ~HPC_flag && PLOT_RMS
        mtit(sprintf('based vibrational component DOI + thermal noise POI %i',iPOI))
        pause(0.1)
    end
    
    [RMS(iPOI).RMSDOIvrpOSCvr,~,RMS(iPOI).Q2DOIvrpOSCvr] = ...
        calcRMS(reconsSignalVRdvpov{iPOI},inputSignal,timeReconS{iPOI},mirror(iPOI),PLOT_RMS);
    if ~HPC_flag && PLOT_RMS
        mtit(sprintf('based vibrational component DOI + OSC vr noise POI %i',iPOI))
        pause(0.1)
    end
    
    [RMS(iPOI).RMSDOIvrpstatnoise,~,RMS(iPOI).Q2DOIvrpstatnoise] = ...
        calcRMS(reconsSignalVRdvpsn{iPOI},inputSignal,timeReconS{iPOI},mirror(iPOI),PLOT_RMS);
    if ~HPC_flag && PLOT_RMS
        mtit(sprintf('based vibrational component DOI + stat noise POI %i',iPOI))
        pause(0.1)
    end
    
    RMS(iPOI).info = num2str(iPOI);
    SNR.info{iPOI} = num2str(iPOI);

end

%% create two extra measurers, one averaged over all POIs and one averaged over POI's highest power DOI 
mVR_POIs = mean(VR.*(-1).^mirror,1);
mVRptherm_POIs = mean(VRpthn.*(-1).^mirror,1);
mVRDOI_POIs = mean(VRDOI.*(-1).^mirror,1);
mVRdv_POIs = mean(VRdv.*(-1).^mirror,1);
mVRdvpthn_POIs = mean(VRdvpthn.*(-1).^mirror,1);
mVRan_POIs = mean(VRan.*(-1).^mirror,1);
mVRanpthn_POIs = mean(VRanpthn.*(-1).^mirror,1);
mVRov_POIs = mean(VRov.*(-1).^mirror,1);
mVRsn_POIs = mean(VRsn.*(-1).^mirror,1);
mVRdvpov_POIs = mean(VRdvpov.*(-1).^mirror,1);
mVRdvpsn_POIs = mean(VRdvpsn.*(-1).^mirror,1);


max_PDOI = max(PowerDOIvr);
idx_PsO = find(PowerDOIvr>max_PDOI/2); % indices of (P)OIs with (s)ame (O)rder of magintude
info_str = ['POIs included:',sprintf(' POI %i,',idx_PsO)]; info_str(end)=[];
if Display
    disp(info_str)
end
mVR_PsO = mean(VR(idx_PsO,:).*(-1).^mirror(idx_PsO),1);
mVRptherm_PsO = mean(VRpthn(idx_PsO,:).*(-1).^mirror(idx_PsO),1);
mVRDOI_PsO = mean(VRDOI(idx_PsO,:).*(-1).^mirror(idx_PsO),1);
mVRdv_PsO = mean(VRdv(idx_PsO,:).*(-1).^mirror(idx_PsO),1);
mVRdvpthn_PsO = mean(VRdvpthn(idx_PsO,:).*(-1).^mirror(idx_PsO),1);
mVRan_PsO = mean(VRan(idx_PsO,:).*(-1).^mirror(idx_PsO),1);
mVRanpthn_PsO = mean(VRanpthn(idx_PsO,:).*(-1).^mirror(idx_PsO),1);
mVRov_PsO = mean(VRov(idx_PsO,:).*(-1).^mirror(idx_PsO),1);
mVRsn_PsO = mean(VRsn(idx_PsO,:).*(-1).^mirror(idx_PsO),1);
mVRdvpov_PsO = mean(VRdvpov(idx_PsO,:).*(-1).^mirror(idx_PsO),1);
mVRdvpsn_PsO = mean(VRdvpsn(idx_PsO,:).*(-1).^mirror(idx_PsO),1);

mVR = [mVR_POIs; mVR_PsO];
mVRpthn = [mVRptherm_POIs;mVRptherm_PsO];
mVRDOI = [mVRDOI_POIs; mVRDOI_PsO];
mVRdv = [mVRdv_POIs; mVRdv_PsO];
mVRdvpthn = [mVRdvpthn_POIs; mVRdvpthn_PsO];
mVRan = [mVRan_POIs; mVRan_PsO];
mVRanpthn = [mVRanpthn_POIs; mVRanpthn_PsO];
mVRov = [mVRov_POIs; mVRov_PsO];
mVRsn = [mVRsn_POIs; mVRsn_PsO];
mVRdvpov = [mVRdvpov_POIs; mVRdvpov_PsO];
mVRdvpsn = [mVRdvpsn_POIs; mVRdvpsn_PsO];



%reconstruct based on complete signal    
[reconsSignalmVR,~,timeReconmS] = calcFourierGetsignal(mVR,Tsim,1:2,fus,fbandwidth,inputSignal,mirror,theta_global,PLOT,'POIstoPLOT',[1,2]);
if ~HPC_flag && PLOT
    mtit('based on all dipoles')
    pause(0.1)
end
%reconstruct based on complete signal plus thermal   
[reconsSignalmVRpthn] = calcFourierGetsignal(mVRpthn,Tsim,1:2,fus,fbandwidth,inputSignal,mirror,theta_global,PLOT,'POIstoPLOT',[1,2]);
if ~HPC_flag && PLOT
    mtit('based on signal all dipoles + thermal noise')
    pause(0.1)
end
% reconstruct based on only DOI
[reconsSignalmVRDOI] = calcFourierGetsignal(mVRDOI,Tsim,1:2,fus,fbandwidth,inputSignal,mirror,theta_global,0);
% reconstruct based on only DOIvibr
[reconsSignalmVRdv,powerSignalmVRdv] = calcFourierGetsignal(mVRdv,Tsim,1:2,fus,fbandwidth,inputSignal,mirror,theta_global,0);
% reconstruct based on only DOIvibrptherm
[reconsSignalmVRdvpthn] = calcFourierGetsignal(mVRdvpthn,Tsim,1:2,fus,fbandwidth,inputSignal,mirror,theta_global,0);
% reconstruct based on all noise
[reconsSignalmVRan,powerSignalmVRan] = calcFourierGetsignal(mVRan,Tsim,1:2,fus,fbandwidth,inputSignal,mirror,theta_global,0);
% reconstruct based on all noise pthn
[reconsSignalmVRanpthn,powerSignalmVRanpthn] = calcFourierGetsignal(mVRanpthn,Tsim,1:2,fus,fbandwidth,inputSignal,mirror,theta_global,0);
% reconstruct based on static noise
[reconsSignalmVRsn,powerSignalmVRsn] = calcFourierGetsignal(mVRsn,Tsim,1:2,fus,fbandwidth,inputSignal,mirror,theta_global,0);
% reconstruct based on OSCvibr
[reconsSignalmVRov,powerSignalmVRov] = calcFourierGetsignal(mVRov,Tsim,1:2,fus,fbandwidth,inputSignal,mirror,theta_global,0);
% reconstruct based on DOI vibr + OSCvibr
[reconsSignalmVRdvpov] = calcFourierGetsignal(mVRdvpov,Tsim,1:2,fus,fbandwidth,inputSignal,mirror,theta_global,0);
% reconstruct based on DOI vibr+ statnoise
[reconsSignalmVRdvpsn] = calcFourierGetsignal(mVRdvpsn,Tsim,1:2,fus,fbandwidth,inputSignal,mirror,theta_global,0);


mVR_strs = {'meanAllPOIs','meanPOIs_sameO'};

for iPOI=1:length(reconsSignalmVR)
    % SNR
    PowerDOIvr(NPOIs+iPOI) = calcPower(reconsSignalmVRdv{iPOI});
    PowerNoiseAll(NPOIs+iPOI) = calcPower(reconsSignalmVRan{iPOI});
    PowerNoiseAllpthn(NPOIs+iPOI) = calcPower(reconsSignalmVRanpthn{iPOI});
    PowerNoiseOSCvr(NPOIs+iPOI) = calcPower(reconsSignalmVRov{iPOI});
    PowerNoiseStatic(NPOIs+iPOI) = calcPower(reconsSignalmVRsn{iPOI});
   
    SNR.DOIvr_NoiseAll(1,NPOIs+iPOI) = 10*log10(PowerDOIvr(NPOIs+iPOI)/PowerNoiseAll(NPOIs+iPOI));
    SNR.DOIvr_NoiseAll(2,NPOIs+iPOI) = 10*log10(powerSignalmVRdv(iPOI)/powerSignalmVRan(iPOI));
    SNR.DOIvr_NoiseAllptherm(1,NPOIs+iPOI) = 10*log10(PowerDOIvr(NPOIs+iPOI)/PowerNoiseAllpthn(NPOIs+iPOI));
    SNR.DOIvr_NoiseAllptherm(2,NPOIs+iPOI) = 10*log10(powerSignalmVRdv(iPOI)/powerSignalmVRanpthn(iPOI));
    SNR.DOIvr_NoiseOSCvr(1,NPOIs+iPOI) = 10*log10(PowerDOIvr(NPOIs+iPOI)/PowerNoiseOSCvr(NPOIs+iPOI));
    SNR.DOIvr_NoiseOSCvr(2,NPOIs+iPOI) = 10*log10(powerSignalmVRdv(iPOI)/powerSignalmVRov(iPOI));
    SNR.DOIvr_NoiseStatic(1,NPOIs+iPOI) = 10*log10(PowerDOIvr(NPOIs+iPOI)/PowerNoiseStatic(NPOIs+iPOI));
    SNR.DOIvr_NoiseStatic(2,NPOIs+iPOI) = 10*log10(powerSignalmVRdv(iPOI)/powerSignalmVRsn(iPOI));
    
    %RMS
    [RMS(NPOIs+iPOI).RMS,~,RMS(NPOIs+iPOI).Q2] = ...
        calcRMS(reconsSignalmVR{iPOI},inputSignal,timeReconmS{iPOI},mirror(1),PLOT_RMS);
    if ~HPC_flag && PLOT_RMS
        mtit(sprintf('based on all dipoles %s',mVR_strs{iPOI}))
        pause(0.1)
    end
    [RMS(NPOIs+iPOI).RMSptherm,~,RMS(NPOIs+iPOI).Q2ptherm] = ...
        calcRMS(reconsSignalmVRpthn{iPOI},inputSignal,timeReconmS{iPOI},mirror(1),PLOT_RMS);
    if ~HPC_flag && PLOT_RMS
        mtit(sprintf('based on all dipoles + thermal noise, %s',mVR_strs{iPOI}))
        pause(0.1)
    end
    
    [RMS(NPOIs+iPOI).RMSDOI,~,RMS(NPOIs+iPOI).Q2DOI] = ...
        calcRMS(reconsSignalmVRDOI{iPOI},inputSignal,timeReconmS{iPOI},mirror(1),PLOT_RMS);
    if ~HPC_flag && PLOT_RMS
        mtit(sprintf('based on DOI, %s',mVR_strs{iPOI}))
        pause(0.1)
    end
    [RMS(NPOIs+iPOI).RMSDOIvr,~,RMS(NPOIs+iPOI).Q2DOIvr] = ...
        calcRMS(reconsSignalmVRdv{iPOI},inputSignal,timeReconmS{iPOI},mirror(1),PLOT_RMS);
    if ~HPC_flag && PLOT_RMS
        mtit(sprintf('based vibrational component DOI, %s',mVR_strs{iPOI}))
        pause(0.1)
    end
    [RMS(NPOIs+iPOI).RMSDOIvrpthn,~,RMS(NPOIs+iPOI).Q2DOIvrpthn] = ...
        calcRMS(reconsSignalmVRdvpthn{iPOI},inputSignal,timeReconmS{iPOI},mirror(1),PLOT_RMS);
    if ~HPC_flag && PLOT_RMS
        mtit(sprintf('based vibrational component DOI + thermal noise, %s',mVR_strs{iPOI}))
        pause(0.1)
    end
    
    [RMS(NPOIs+iPOI).RMSDOIvrpOSCvr,~,RMS(NPOIs+iPOI).Q2DOIvrpOSCvr] = ...
        calcRMS(reconsSignalmVRdvpov{iPOI},inputSignal,timeReconmS{iPOI},mirror(1),PLOT_RMS);
    if ~HPC_flag && PLOT_RMS
        mtit(sprintf('based vibrational component DOI + OSC vibr, %s',mVR_strs{iPOI}))
        pause(0.1)
    end
    
    [RMS(NPOIs+iPOI).RMSDOIvrpstatnoise,~,RMS(NPOIs+iPOI).Q2DOIvrpstatnoise] = ...
        calcRMS(reconsSignalmVRdvpsn{iPOI},inputSignal,timeReconmS{iPOI},mirror(1),PLOT_RMS);
    if ~HPC_flag && PLOT_RMS
        mtit(sprintf('based vibrational component DOI + static noise, %s',mVR_strs{iPOI}))
        pause(0.1)
    end
    
    RMS(NPOIs+iPOI).info = mVR_strs{iPOI};
    SNR.info{NPOIs+iPOI} = mVR_strs{iPOI};

end


timeReconS_f = horzcat(timeReconS,timeReconmS);
reconsSignalVR_f = horzcat(reconsSignalVR,reconsSignalmVR);
reconsSignalVRDOI_f = horzcat(reconsSignalVRDOI,reconsSignalmVRDOI);
reconsSignalVRdv_f = horzcat(reconsSignalVRdv,reconsSignalmVRdv);
reconsSignalVRpthn_f = horzcat(reconsSignalVRpthn,reconsSignalmVRpthn);
reconsSignalVRdvpov_f = horzcat(reconsSignalVRdvpov,reconsSignalmVRdvpov);
reconsSignalVRdvpsn_f = horzcat(reconsSignalVRdvpsn,reconsSignalmVRdvpsn);
VR_f = vertcat(VR,mVR);
VRpthn_f = vertcat(VRpthn,mVRpthn);
VRDOI_f = vertcat(VRDOI,mVRDOI);
VRDOIvibr_f = vertcat(VRdv,mVRdv);



    function Power = calcPower(Signal)
        s0dft = fft(Signal);
        psds0 = (1/(length(Signal))) * abs(s0dft).^2;
        Power = sum(psds0);
    end

    


end