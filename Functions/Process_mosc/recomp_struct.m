function [RMSall,SNRall,Outall,TSTOPall] = recomp_struct(RMS_mat,RMSptherm_mat,RMSDOI_mat,RMSDOIvr_mat,RMSDOIvrptherm_mat,...
    RMSDOIvrpOSCvr_mat,RMSDOIvrpstatnoise_mat,RMS_info,...
    Q2_mat,Q2ptherm_mat,Q2DOI_mat,Q2DOIvr_mat,Q2DOIvrptherm_mat,Q2DOIvrpOSCvr_mat,Q2DOIvrpstatnoise_mat,...
    rSVR_mat,rSVRDOIvr_mat,rSVRDOIptherm_mat,rSVRdvpov_mat,rSVRdvpsn_mat,rSVRDOI_mat,trS_mat,...
    SNR_DOIvr_noiseall_mat,SNR_DOIvr_noiseallptherm_mat,SNR_DOIvr_noiseoscvr_mat,SNR_DOIvr_noisestat_mat,...
    SNR_info,VR_mat,VRptherm_mat,VRDOI_mat,VRDOIstat_mat,VROSC_mat,VROSCstat_mat,VRstatnoise_mat,StartSim_mat,Itgen_mat,...
    VRgen_mat,endsim_mat,signalreconstruction_mat,POIs,Input,Outall_in,nos4l_flag)

% for memory reasons, at the end of calcErrorGrid, the structure
% is decomposed into arrays => lower memory to save. Recompose back to
% structure in order to analyse the results


liPOIs = size(POIs,1);
[liAmp,lidps,liPOIs_RMS] = size(RMS_mat);

for iAmp = 1:liAmp
    for idps = 1:lidps
        for iPOIs = 1:liPOIs_RMS
            
            RMSall(iAmp,idps).RMS(iPOIs).RMS = RMS_mat(iAmp,idps,iPOIs);
            RMSall(iAmp,idps).RMS(iPOIs).Q2 = Q2_mat(iAmp,idps,iPOIs);
            RMSall(iAmp,idps).RMS(iPOIs).RMSptherm = RMSptherm_mat(iAmp,idps,iPOIs);
            RMSall(iAmp,idps).RMS(iPOIs).Q2ptherm = Q2ptherm_mat(iAmp,idps,iPOIs);
            RMSall(iAmp,idps).RMS(iPOIs).RMSDOI = RMSDOI_mat(iAmp,idps,iPOIs);            
            RMSall(iAmp,idps).RMS(iPOIs).Q2DOI = Q2DOI_mat(iAmp,idps,iPOIs);
            RMSall(iAmp,idps).RMS(iPOIs).RMSDOIvr = RMSDOIvr_mat(iAmp,idps,iPOIs);            
            RMSall(iAmp,idps).RMS(iPOIs).Q2DOIvr = Q2DOIvr_mat(iAmp,idps,iPOIs);            
            RMSall(iAmp,idps).RMS(iPOIs).RMSDOIvrpthn = RMSDOIvrptherm_mat(iAmp,idps,iPOIs);
            RMSall(iAmp,idps).RMS(iPOIs).Q2DOIvrpthn = Q2DOIvrptherm_mat(iAmp,idps,iPOIs);
            RMSall(iAmp,idps).RMS(iPOIs).RMSDOIvrpOSCvr = RMSDOIvrpOSCvr_mat(iAmp,idps,iPOIs);
            RMSall(iAmp,idps).RMS(iPOIs).Q2DOIvrpOSCvr = Q2DOIvrpOSCvr_mat(iAmp,idps,iPOIs);            
            RMSall(iAmp,idps).RMS(iPOIs).RMSDOIvrpstatnoise = RMSDOIvrpstatnoise_mat(iAmp,idps,iPOIs);
            RMSall(iAmp,idps).RMS(iPOIs).Q2DOIvrpstatnoise = Q2DOIvrpstatnoise_mat(iAmp,idps,iPOIs);            
            RMSall(iAmp,idps).RMS(iPOIs).info = RMS_info{iAmp,idps,iPOIs};
            
            %recompose Out reconsSignal    
            Outall(iAmp,idps).Out.reconsSignalVR(iPOIs) = {squeeze(rSVR_mat(iAmp,idps,iPOIs,:))'};
            Outall(iAmp,idps).Out.reconsSignalVRDOI(iPOIs) = {squeeze(rSVRDOI_mat(iAmp,idps,iPOIs,:))'};
            Outall(iAmp,idps).Out.reconsSignalVRDOIvibr(iPOIs) = {squeeze(rSVRDOIvr_mat(iAmp,idps,iPOIs,:))'};
            Outall(iAmp,idps).Out.reconsSignalVRDOIptherm(iPOIs) = {squeeze(rSVRDOIptherm_mat(iAmp,idps,iPOIs,:))'};
            Outall(iAmp,idps).Out.reconsSignalVRdvpov(iPOIs) = {squeeze(rSVRdvpov_mat(iAmp,idps,iPOIs,:))'};
            Outall(iAmp,idps).Out.reconsSignalVRdvpsn(iPOIs) = {squeeze(rSVRdvpsn_mat(iAmp,idps,iPOIs,:))'};
            Outall(iAmp,idps).Out.timereconS(iPOIs) = {squeeze(trS_mat(iAmp,idps,iPOIs,:))'};
            
            
            % Decompose SNR structure
            SNRall(iAmp,idps).SNR.DOIvr_NoiseAll = squeeze(SNR_DOIvr_noiseall_mat(iAmp,idps,:,:));
            SNRall(iAmp,idps).SNR.DOIvr_NoiseAllptherm = squeeze(SNR_DOIvr_noiseallptherm_mat(iAmp,idps,:,:));
            SNRall(iAmp,idps).SNR.DOIvr_NoiseOSCvr = squeeze(SNR_DOIvr_noiseoscvr_mat(iAmp,idps,:,:));
            SNRall(iAmp,idps).SNR.DOIvr_NoiseStatic = squeeze(SNR_DOIvr_noisestat_mat(iAmp,idps,:,:));
            SNRall(iAmp,idps).SNR.info(:) = squeeze(SNR_info(iAmp,idps,:));
            
            % Recompose Out
            
            Outall(iAmp,idps).Out.VR = squeeze(VR_mat(iAmp,idps,:,:));
            Outall(iAmp,idps).Out.VRDOI = squeeze(VRDOI_mat(iAmp,idps,:,:));
            Outall(iAmp,idps).Out.VRptherm = squeeze(VRptherm_mat(iAmp,idps,:,:));            
            Outall(iAmp,idps).Out.VRDOIstat = squeeze(VRDOIstat_mat(iAmp,idps,:,:));
            Outall(iAmp,idps).Out.VROSC = squeeze(VROSC_mat(iAmp,idps,:,:));
            Outall(iAmp,idps).Out.VROSCstat = squeeze(VROSCstat_mat(iAmp,idps,:,:));
            Outall(iAmp,idps).Out.VRstatnoise = squeeze(VRstatnoise_mat(iAmp,idps,:,:));
            
            %Recompose TStOP
            
            TSTOPall(iAmp,idps).TSTOP.startSim = squeeze(StartSim_mat(iAmp,idps,:))';
            TSTOPall(iAmp,idps).TSTOP.Itgen = squeeze(Itgen_mat(iAmp,idps,:))' ;
            TSTOPall(iAmp,idps).TSTOP.VRgen = squeeze(VRgen_mat(iAmp,idps,:))' ;
            TSTOPall(iAmp,idps).TSTOP.endsim = squeeze(endsim_mat(iAmp,idps,:))';
            TSTOPall(iAmp,idps).TSTOP.signalreconstruction = squeeze(signalreconstruction_mat(iAmp,idps,:))';
            
            if ~isempty(Outall_in(iAmp,idps).Out)
                Outall(iAmp,idps).Out.CSource = Outall_in(iAmp,idps).Out.CSource;
                Outall(iAmp,idps).Out.CSink = Outall_in(iAmp,idps).Out.CSink;
                Outall(iAmp,idps).Out.POIs = Outall_in(iAmp,idps).Out.POIs;
                Outall(iAmp,idps).Out.OSCindices = Outall_in(iAmp,idps).Out.OSCindices;
                Outall(iAmp,idps).Out.randomseed = Outall_in(iAmp,idps).Out.randomseed;
                Outall(iAmp,idps).Out.extra_opt = Outall_in(iAmp,idps).Out.extra_opt;
                Outall(iAmp,idps).Out.Param = Outall_in(iAmp,idps).Out.Param;
            end
            
            
            if ~nos4l_flag
                if ~isempty(Outall_in(iAmp,idps).Out)
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

end