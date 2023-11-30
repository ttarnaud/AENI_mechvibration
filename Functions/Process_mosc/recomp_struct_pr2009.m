function [RMSall,SNRall,Outall,TSTOPall] = recomp_struct_pr2009(RMS_mat,RMSDOI_mat,Q2_mat,Q2DOI_mat,RMSptherm_mat,...
    RMSDOIptherm_mat,Q2ptherm_mat,Q2DOIptherm_mat,RMS_info,rSVR_mat,rSVRDOI_mat,trS_mat,...
    SNR_DOI_noiseall_mat,SNR_DOI_noiseallptherm_mat,SNR_DOI_noiseosc_mat,SNR_DOI_noisestat_mat,...
    SNR_info,VR_mat,VRptherm_mat,VRDOI_mat,VRnoise_stat_mat,StartSim_mat,Itgen_mat,...
    VRgen_mat,endsim_mat,signalreconstruction_mat,POIs,Input,Outall_in,nos4l_flag)

% for memory reasons, at the end of invest biological noise, the structure
% is decomposed into arrays => lower memory to save. Recompose back to
% structure in order to analyse the results

% use this recomposition method for results prior september 2020; 
% the averaged results (over all POIs and pois with same oreder of
% magnitude stored in diffrent way)

liPOIs = size(POIs,1);
[liAmp,lidps,liPOIs_RMS] = size(RMS_mat);

for iAmp = 1:liAmp
    for idps = 1:lidps
        for iPOIs = 1:liPOIs_RMS
            
            RMSall(iAmp,idps).RMS(iPOIs).RMS = RMS_mat(iAmp,idps,iPOIs);
            RMSall(iAmp,idps).RMS(iPOIs).RMSDOI = RMSDOI_mat(iAmp,idps,iPOIs);
            RMSall(iAmp,idps).RMS(iPOIs).Q2 = Q2_mat(iAmp,idps,iPOIs);
            RMSall(iAmp,idps).RMS(iPOIs).Q2DOI = Q2DOI_mat(iAmp,idps,iPOIs);
            RMSall(iAmp,idps).RMS(iPOIs).RMSptherm = RMSptherm_mat(iAmp,idps,iPOIs);
            RMSall(iAmp,idps).RMS(iPOIs).RMSDOIptherm = RMSDOIptherm_mat(iAmp,idps,iPOIs);
            RMSall(iAmp,idps).RMS(iPOIs).Q2ptherm = Q2ptherm_mat(iAmp,idps,iPOIs);
            RMSall(iAmp,idps).RMS(iPOIs).Q2DOIptherm = Q2DOIptherm_mat(iAmp,idps,iPOIs);
            RMSall(iAmp,idps).RMS(iPOIs).info = RMS_info{iAmp,idps,iPOIs};
            if iPOIs<=liPOIs
                %Decompose Out reconsSignal
                Outall(iAmp,idps).Out.reconsSignalVR(iPOIs) = {squeeze(rSVR_mat(iAmp,idps,iPOIs,:))'};
                Outall(iAmp,idps).Out.reconsSignalVRDOI(iPOIs) = {squeeze(rSVRDOI_mat(iAmp,idps,iPOIs,:))'};
                Outall(iAmp,idps).Out.timereconS(iPOIs) = {squeeze(trS_mat(iAmp,idps,iPOIs,:))'};
            end
            
            % Decompose SNR structure
            SNRall(iAmp,idps).SNR.DOI_noiseall = squeeze(SNR_DOI_noiseall_mat(iAmp,idps,:,:));
            SNRall(iAmp,idps).SNR.DOI_noiseallptherm = squeeze(SNR_DOI_noiseallptherm_mat(iAmp,idps,:,:));
            SNRall(iAmp,idps).SNR.DOI_noiseosc = squeeze(SNR_DOI_noiseosc_mat(iAmp,idps,:,:));
            SNRall(iAmp,idps).SNR.DOI_noisestat = squeeze(SNR_DOI_noisestat_mat(iAmp,idps,:,:));
            SNRall(iAmp,idps).SNR.info(:) = squeeze(SNR_info(iAmp,idps,:));
            
            % Decompose Out
            Outall(iAmp,idps).Out.VR = squeeze(VR_mat(iAmp,idps,1:liPOIs,:));
            Outall(iAmp,idps).Out.VRptherm = squeeze(VRptherm_mat(iAmp,idps,1:liPOIs,:));
            Outall(iAmp,idps).Out.VRDOI = squeeze(VRDOI_mat(iAmp,idps,1:liPOIs,:));
            Outall(iAmp,idps).Out.VRnoise_stat = squeeze(VRnoise_stat_mat(iAmp,idps,1:liPOIs,:));
            Outall(iAmp,idps).Out.mVR_POI = squeeze(VR_mat(iAmp,idps,liPOIs+1,:))';
            Outall(iAmp,idps).Out.mVR_PsO = squeeze(VR_mat(iAmp,idps,liPOIs+2,:))';
            Outall(iAmp,idps).Out.mVRDOI_POI = squeeze(VRDOI_mat(iAmp,idps,liPOIs+1,:))' ;
            Outall(iAmp,idps).Out.mVRDOI_PsO = squeeze(VRDOI_mat(iAmp,idps,liPOIs+2,:))';
            Outall(iAmp,idps).Out.mVRptherm_POI = squeeze(VRptherm_mat(iAmp,idps,liPOIs+1,:))';
            Outall(iAmp,idps).Out.mVRptherm_PsO = squeeze(VRptherm_mat(iAmp,idps,liPOIs+2,:))';
            Outall(iAmp,idps).Out.mVRnoise_stat_POIs = squeeze(VRnoise_stat_mat(iAmp,idps,liPOIs+1,:))';
            Outall(iAmp,idps).Out.mVRnoise_stat_PsO = squeeze(VRnoise_stat_mat(iAmp,idps,liPOIs+2,:))';
            
            
            
            %Decompose Out reconsSignal
            Outall(iAmp,idps).Out.reconsSignalmVR_POI(1) = {squeeze(rSVR_mat(iAmp,idps,liPOIs+1,:))'};
            Outall(iAmp,idps).Out.reconsSignalmVR_PsO(1) = {squeeze(rSVR_mat(iAmp,idps,liPOIs+2,:))'};
            Outall(iAmp,idps).Out.reconsSignalmVRDOI_POI(1) = {squeeze(rSVRDOI_mat(iAmp,idps,liPOIs+1,:))'};
            Outall(iAmp,idps).Out.reconsSignalmVRDOI_PsO(1) = {squeeze(rSVRDOI_mat(iAmp,idps,liPOIs+2,:))'};
            Outall(iAmp,idps).Out.timereconS_POI(1) = {squeeze(trS_mat(iAmp,idps,liPOIs+1,:))'};
            Outall(iAmp,idps).Out.timereconS_PsO(1) = {squeeze(trS_mat(iAmp,idps,liPOIs+2,:))'};
            
            %Decompose TStOP
            
            TSTOPall(iAmp,idps).TSTOP.startSim = squeeze(StartSim_mat(iAmp,idps,:))';
            TSTOPall(iAmp,idps).TSTOP.Itgen = squeeze(Itgen_mat(iAmp,idps,:))' ;
            TSTOPall(iAmp,idps).TSTOP.VRgen = squeeze(VRgen_mat(iAmp,idps,:))' ;
            TSTOPall(iAmp,idps).TSTOP.endsim = squeeze(endsim_mat(iAmp,idps,:))';
            TSTOPall(iAmp,idps).TSTOP.signalreconstruction = squeeze(signalreconstruction_mat(iAmp,idps,:))';
            
            Outall(iAmp,idps).Out.CSource = Outall_in(iAmp,idps).Out.CSource;
            Outall(iAmp,idps).Out.CSink = Outall_in(iAmp,idps).Out.CSink;
            Outall(iAmp,idps).Out.POIs = Outall_in(iAmp,idps).Out.POIs;
            Outall(iAmp,idps).Out.OSCindices = Outall_in(iAmp,idps).Out.OSCindices;
            
            if ~nos4l_flag
            Outall(iAmp,idps).Out.vampx = Outall_in(iAmp,idps).Out.vampx;
            Outall(iAmp,idps).Out.vampy = Outall_in(iAmp,idps).Out.vampy;
            Outall(iAmp,idps).Out.vampz = Outall_in(iAmp,idps).Out.vampz;
            Outall(iAmp,idps).Out.phasex = Outall_in(iAmp,idps).Out.phasex;
            Outall(iAmp,idps).Out.vphasey = Outall_in(iAmp,idps).Out.vphasey;
            Outall(iAmp,idps).Out.vphasez = Outall_in(iAmp,idps).Out.vphasex;
                
            end
            
        end
    end
end

end