function [requiredAmp,requiredAmplower,itermat,cSC_mode] = ...
    cSC_cEG(Input,Amp,kAus,kAus_con,vibrField,flaglim,...
    evalMode,maxTime,maxTime_overall,timeGrid,...
    Totaldps,dps_run,RMS_mat,Threshold,maxIter,input_set_init,...
    SettingsStr,Input_str,s)
% calculate strength curve (cSC) =  minimal amplitude our spatial constant
% in order to get a RMSE lower than a given percentage. This function is
% used in the calcErrorGrid master (cEG)

requiredAmp = [];
requiredAmplower = [];
itermat = []
cSC_mode = 'amp';
if ~any(Amp-Amp(1))  && strcmpi(vibrField,'art')
    if length(kAus)>1
    cSC_mode = 'kAus';
    else
        warning('calc strength mode is kAus but only one value?')
    end
end
if isfield(Input,'cSC_mode')
    cSC_mode = Input.cSC_mode;
end

switch lower(cSC_mode)
    case 'amp'
        optim_val = Amp;
        flagfun = @(a,b) abs(a-b);
        if ~isempty(kAus) 
            if ~(length(kAus)==length(Totaldps)||length(kAus)==1)
            error('length kAus must be same as length Totaldps or equal to one in cSC_mode amp')
            end
        end
    case 'kaus'
        optim_val = kAus;
        flaglim = 0.0001;
        flagfun = @(a,b) abs(1./a-1./b);
        if length(kAus)~=length(Amp)
            error('length Amp and kAus must be same in cSC_mode kaus')
        end
    otherwise
        error('false cSC_mode')
end
fprintf('\n start calculation of strength Curve\n');
if strcmpi(evalMode,'hpc')
    maxTime = 0.9*(maxTime_overall-timeGrid);
end

for idps = 1:length(Totaldps)
    %time limit!
    maxTime_dprun = floor(maxTime.*Totaldps(idps)./sum(Totaldps));
    %reset flags
    iter = 0;
    findHigherAmpiter = 0;
    startRun = tic; %start timer
    
    RMS_mat_sel = squeeze(RMS_mat(:,idps,:));
    RMS_minPOI = min(RMS_mat_sel,[],2);
    superval = optim_val(RMS_minPOI<=Threshold); %amplitude that causes total error to be lower than thershold
    
    if isempty(superval)
        superval = max(optim_val);
        newval = superval;
        findhigheramp = 1;
    else
        superval = superval(1);
        findhigheramp = 0;
    end
    
    lowerval = optim_val(RMS_minPOI>=Threshold); %amplitude that causes total error to be higher than thershold
    if isempty(lowerval)
        lowerval = superval(1);
    else
        lowerval = lowerval(end);
    end
    flag = findhigheramp || flagfun(superval,lowerval)>flaglim;
    iterTimeStop = zeros(1,maxIter);
    while flag && iter<maxIter && findHigherAmpiter<10
        iterTimeStart = tic;
        % reset random generator back to original seed => threshold search not
        % dependend on random generation each iteration
        rng(s)
        if findhigheramp
            newval = newval*10;
            iter = 0;
            findHigherAmpiter = findHigherAmpiter +1;
        else
            newval = 10^((log10(superval)+log10(lowerval))/2);
        end
        fprintf('\n%s\n %s: %5.2e, idps: %i/%i, iter: %i\n',datestr(now),cSC_mode,newval,idps,length(Totaldps),iter)
        
        switch lower(vibrField)
            case 's4l'
                Amp_input = newval;
                Input_AdaptVamp.Aus = newval;
                k_Aus = [];
            case 'art'
                switch lower(cSC_mode)
                    case 'amp'
                        
                        switch lower(kAus_con)
                            case {'none','amp'}
                                k_Aus = kAus(1);
                                if length(kAus)>1 && any(kAus-kAus(1))
                                    warning('Amp optimisation with mulitple kAus?')
                                end
                            case 'totaldps'
                                k_Aus = kAus(idps);
                            otherwise
                                error('false kAus_con input')
                        end                                              
                        Amp_input = newval;
                    case 'kaus'
                        k_Aus = newval;
                        Amp_input = Amp(1);
                    otherwise
                        error('false cSC_mode')
                end
                
            case 'none'
                k_Aus = [];
                Amp_input = newval;
        end
        input_set = horzcat({'Totaldps',Totaldps(idps),'dps_run',dps_run(idps),'Aus',Amp_input,...
                'k_Aus',k_Aus},input_set_init);
 
        [RMS,SNR,TSTOP,Out] = investBiologicalNoise(input_set{:});
        RMSn = min([RMS(:).RMS]);
        if RMSn<=Threshold
            superval = newval;
            findhigheramp = 0;
        else
            lowerval = newval;
        end
        timethisdps = toc(startRun);
        iterTimeStop(iter+1) = toc(iterTimeStart);
        flag = flagfun(superval,lowerval)>flaglim & (timethisdps+max(iterTimeStop))<maxTime_dprun;
        iter = iter+1;
        
        
    end
    
    if flagfun(superval,lowerval)<=flaglim
        fprintf('\n end iteration due to reaching flaglimit: %5.2e\n',flaglim)
    elseif(timethisdps+max(iterTimeStop))>=maxTime_dprun
        fprintf('\n end iteration due to reaching timelimit: %5.2f\n',maxTime_dprun)
    elseif iter>=maxIter
        fprintf('\n end iteration due to reaching iteration limit: %i\n',maxIter)
    else
        disp('end iteration for some reason')
    end
    
    requiredAmp(idps) = superval;
    requiredAmplower(idps) = lowerval;
    itermat(idps) = iter;
    
    %save intermediate results
    idxdoublepoint = find(SettingsStr==':');
    if isempty(idxdoublepoint)
        save_str_int = SettingsStr;
    else
        save_str_int = SettingsStr([1:idxdoublepoint-1,idxdoublepoint+2:end]);
    end
    % %         savename = sprintf('Results_StrengthCurve%s_%s_%s_%s.mat',savename_add,save_str_int,...
    % %             Input_str(regexp(Input_str,'_sM')+3:end),datestr(now,'mm-dd-yy_HH'));
    savename = sprintf('Results_StrengthCurve_%s_%s_%s.mat',save_str_int,...
        Input_str(regexp(Input_str,'_cEG')+5:end-2),datestr(now,'mm-dd-yy_HH'));
    save(savename,'requiredAmp','requiredAmplower','itermat','cSC_mode');
    
end


end