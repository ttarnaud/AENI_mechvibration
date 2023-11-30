function [CSource,CSink,OSCindices,DOIindice,USwave_all,Dirus_fun,Phaseus_fun,Aus_fun,...
    Aus,Phaseus,Dirus,wus,USwave_allsame_flag,varstoclear,plotUSwaves_flag,vampx_out,vampy_out,vampz_out] =...
    declare_vibrationfunc(posDp,d,CSource,CSink,OrienDipole,Dirus,Aus,wus,Phaseus,Dirus_way,Aus_way,k_Aus,...
    USwave_varlengths,eval_method,all_osc_flag,S4l_flag,S4l_Input,interp3_method,AdaptVamp_flag,Aus_AVa,vmax_AVa,...
    sel_oscil_from_gendPos_flag,ROI_OSC,plotUSwaves_flag,coll_singlePeriod_flag,rel_singlePeriod_flag,Input_VRsP,Display)
% in this function we determine the oscilating dipoles and assign a
% corresponing wave function US_wave which alters te position of the dipole
% in function of time. ==> vibrates


varstoclear={};
vampx_out = [];
vampy_out = [];
vampz_out = [];

%% get OSCindices
if sel_oscil_from_gendPos_flag
    
    % only assign vibration to dipoles in certain region specified by
    % ROI_OSC
    CSource = vertcat(posDp+d/2*OrienDipole,CSource);
    CSink = vertcat(posDp-d/2*OrienDipole,CSink);
    OSCindices = find(inhull(CSource,ROI_OSC)&inhull(CSink,ROI_OSC));
elseif S4l_flag
    
    % Assign vibration based on pressure profile obtained in S4L
    CSource = vertcat(posDp+d/2.*OrienDipole,CSource);
    CSink = vertcat(posDp-d/2.*OrienDipole,CSink);
    posDPall = (CSource+CSink)/2;
    
    % check if all dipoles are in brain simulation domain
    if any(vecnorm(CSource,2,2)>S4l_Input.RBrain) || any(vecnorm(CSink,2,2)>S4l_Input.RBrain)
        error('dipoles are not in simulation domain s4l')
    end
    
    % Unit conversion
    switch S4l_Input.loc.unit
        case 'm'
            f_loc = 1;
        case 'mm'
            f_loc = 1e-3;
        otherwise
            error('unit not implemented')
    end
    
    S4l_Input.loc.X = S4l_Input.loc.X*f_loc;
    S4l_Input.loc.Y = S4l_Input.loc.Y*f_loc;
    S4l_Input.loc.Z = S4l_Input.loc.Z*f_loc;
    S4l_maxX = max(S4l_Input.loc.X(:)); S4l_minX = min(S4l_Input.loc.X(:));
    S4l_maxY = max(S4l_Input.loc.Y(:)); S4l_minY = min(S4l_Input.loc.Y(:));
    S4l_maxZ = max(S4l_Input.loc.Z(:)); S4l_minZ = min(S4l_Input.loc.Z(:));
    S4l_maxloc = [S4l_maxX, S4l_maxY, S4l_maxZ];
    S4l_minloc = [S4l_minX, S4l_minY, S4l_minZ];
    
    % output of s4l can be cropped for memory reasons find which are
    % contained within boundary
    OSCindices = find(~(any([S4l_maxloc-(CSource+CSink)/2]<0,2)|...
        any([S4l_minloc-(CSource+CSink)/2]>0,2)));
    
    % assign vibration amplitudes and phases to dipoles
    vampx = interp3(S4l_Input.loc.X,S4l_Input.loc.Y,S4l_Input.loc.Z,S4l_Input.veloc.ampx,...
        posDPall(:,1),posDPall(:,2),posDPall(:,3),interp3_method);
    vampy = interp3(S4l_Input.loc.X,S4l_Input.loc.Y,S4l_Input.loc.Z,S4l_Input.veloc.ampy,...
        posDPall(:,1),posDPall(:,2),posDPall(:,3),interp3_method);
    vampz = interp3(S4l_Input.loc.X,S4l_Input.loc.Y,S4l_Input.loc.Z,S4l_Input.veloc.ampz,...
        posDPall(:,1),posDPall(:,2),posDPall(:,3),interp3_method); 
    vphasex = interp3(S4l_Input.loc.X,S4l_Input.loc.Y,S4l_Input.loc.Z,S4l_Input.veloc.phasex,...
        posDPall(:,1),posDPall(:,2),posDPall(:,3),interp3_method);
    vphasey = interp3(S4l_Input.loc.X,S4l_Input.loc.Y,S4l_Input.loc.Z,S4l_Input.veloc.phasey,...
        posDPall(:,1),posDPall(:,2),posDPall(:,3),interp3_method);
    vphasez = interp3(S4l_Input.loc.X,S4l_Input.loc.Y,S4l_Input.loc.Z,S4l_Input.veloc.phasez,...
        posDPall(:,1),posDPall(:,2),posDPall(:,3),interp3_method);
    if coll_singlePeriod_flag
        warning('scaling to SimDipoleOSc not added... how combine to one parameter Aus?, taking norm of vamp values?')
        % should be something of storing Parmaeters above with
        % coll_singlePeriod_flag in Input_VRsP and scaling them when
        % rel_singlePeriod_flag is on Butmemory cost?
        vampx_out = vampx;
        vampy_out = vampy;
        vampz_out = vampz;
    end
    if rel_singlePeriod_flag
        error('scaling to SimDipoleOSc not added... how combine to one parameter Aus?')
        vampx = vampx./Input_VRsP.vampx;
        vampy = vampy./Input_VRsP.vampy;
        vampz = vampz./Input_VRsP.vampz;
    end
    
    if any(strcmpi(eval_method,{'database','database_f','validate'}))
        maxvx = interp3(S4l_Input.loc.X,S4l_Input.loc.Y,S4l_Input.loc.Z,S4l_Input.veloc.maxvx,...
            posDPall(:,1),posDPall(:,2),posDPall(:,3),interp3_method);
        maxvy = interp3(S4l_Input.loc.X,S4l_Input.loc.Y,S4l_Input.loc.Z,S4l_Input.veloc.maxvy,...
            posDPall(:,1),posDPall(:,2),posDPall(:,3),interp3_method);
        maxvz = interp3(S4l_Input.loc.X,S4l_Input.loc.Y,S4l_Input.loc.Z,S4l_Input.veloc.maxvz,...
            posDPall(:,1),posDPall(:,2),posDPall(:,3),interp3_method);
        phasemaxdir = interp3(S4l_Input.loc.X,S4l_Input.loc.Y,S4l_Input.loc.Z,S4l_Input.veloc.phasemaxdir,...
            posDPall(:,1),posDPall(:,2),posDPall(:,3),interp3_method);        
    end
    clear('S4l_Input') %memory reasons
    varstoclear = horzcat(varstoclear,{'S4l_Input'});
else
    % Oscillators are the ones given in posDp
    posOsc = [posDp]; % position of oscilators
    CSource = vertcat(posDp+d/2.*OrienDipole,CSource);
    CSink = vertcat(posDp-d/2.*OrienDipole,CSink);
    OSCindices = 1:size(posOsc,1);
    if all_osc_flag
        OSCindices = 1:size(CSource,1);
    end
end
%%
% Dipole of interest is always the first in the dipole position array
DOIindice = 1;
if Display
    fprintf('\n Dipole of interest is in position %i \n',DOIindice)
end

if ~any(DOIindice==OSCindices)
    warning('DOI not an oscillator')
end
%% Create USwave_all function. This is called for assigning vibration
%direction an amplitude
USwave_allsame_flag = 0;
Dirus_fun = [];
Phaseus_fun = [];
Aus_fun = [];
vibrParam_flag = any(strcmpi(eval_method,{'database','database_f','validate'}));
if ~S4l_flag
    switch lower(Dirus_way)
        case 'dipole_dir'
            if Display
                fprintf('\n Dirus changed from input:  %s \n',Dirus_way)
            end
            Dirus_DOI = Dirus(DOIindice,:);                  %store setting DOI => not altered
            Dirus = CSource(OSCindices,:)-CSink(OSCindices,:);
            Dirus = Dirus./vecnorm(Dirus);
            Dirus(DOIindice,:) = Dirus_DOI;
            USwave_varlengths = [size(Dirus,1), size(Aus,1), length(wus), size(Phaseus,1)];
        case 'radial'
            if Display
                fprintf('\n Dirus changed from input:  %s \n',Dirus_way)
            end
            Dirus_DOI = Dirus(DOIindice,:);
            Dirus = (CSource(OSCindices,:)+CSink(OSCindices,:))/2;
            Dirus = Dirus./vecnorm(Dirus);
            Dirus(DOIindice,:) = Dirus_DOI;
            USwave_varlengths = [size(Dirus,1), size(Aus,1), length(wus), size(Phaseus,1)];
        case 'input'
            if Display
                fprintf('\n Dirus not changed from input')
            end
        otherwise 
            error('false input Dirus_way')

    end
    switch lower(Aus_way)
        case 'kexp'
            % exponential decreasing centered around DOI
            if Display
                fprintf('\n Aus changed from input: %s, k=%3.2e',Aus_way,k_Aus)
            end
            if isrow(Aus)
                Aus = Aus';
            end
            Aus_DOI = Aus(DOIindice,:);
            kfun = @(x)  Aus_DOI*exp(-k_Aus*x);
            Rdiffdppos = vecnorm((CSource(OSCindices,:)+CSink(OSCindices,:)-CSource(DOIindice,:)-CSink(DOIindice,:))/2,2,2);
            Aus = kfun(Rdiffdppos);
            Aus(DOIindice,:) = Aus_DOI;            
            USwave_varlengths = [size(Dirus,1), size(Aus,1), length(wus), size(Phaseus,1)];
        case 'kscale'
            %scaled to value of DOI
            if Display
                fprintf('\n Aus changed from input: %s, k=%3.2e',Aus_way,k_Aus)
            end
            if isrow(Aus)
                Aus = Aus';
            end
            Aus_DOI = Aus(DOIindice,:);
            Aus = Aus_DOI*k_Aus*ones(size(CSource(OSCindices,:),1));
            Aus(DOIindice,:) = Aus_DOI;
            USwave_varlengths = [size(Dirus,1), size(Aus,1), length(wus), size(Phaseus,1)];
        case 'normal'
            if Display
                fprintf('\n Aus not changed from input\n')
            end
        otherwise
            error('false input Aus_way')
    end
    if Display
        disp(' ')
    end
    if rel_singlePeriod_flag
        % scaling Amplitude to strength used in collection phase
        Aus_relsP = Aus./Input_VRsP.Aus;
    end
    % vibration profile not given through S4L
    if length(OSCindices) == 1
        
        % if dimensions of USwaves bigger display warning
        if any(USwave_varlengths-1)
            warning('USwaves more dimension than selected oscillators')
        end
        % because there is only one oscillator we take always first
        % setting set as USwave
        USwave_all = @(t,idx) Dirus(idx==OSCindices,:).*Aus(idx==OSCindices,:).*sin(wus(idx==OSCindices).*t+Phaseus(idx==OSCindices,:));
        USwave_allsame_flag = 1;
        if vibrParam_flag
            Dirus_fun = @(idx) Dirus(idx==OSCindices,:);
            Phaseus_fun = @(idx) Phaseus(idx==OSCindices,:);
            Aus_fun = @(idx) Aus(intersect(idx,OSCindices),:);
        end
        if rel_singlePeriod_flag
            Aus_fun = @(idx) Aus_relsP(intersect(idx,OSCindices),:);
        end
            
    
    elseif isempty(OSCindices)
        
        warning('no oscillator selected')
        USwave_all = @(t,idx) 0*Dirus.*Aus.*sin(wus.*t+Phaseus);
        % eventhough no oscillator selected still declare USwave_all
        plotUSwaves_flag = 0;
        
    else
        % check if dimensions USwave == 1
        if ~any(USwave_varlengths-1)
            warning('all oscillators have same USwave')
            USwave_all = @(t,idx) Dirus.*Aus.*sin(wus.*t+Phaseus);
            USwave_allsame_flag = 1;
        if vibrParam_flag
            Dirus_fun = @(idx) Dirus;
            Phaseus_fun = @(idx) Phaseus;
            Aus_fun = @(idx) Aus.*ones(length(idx),1);
        end
        if rel_singlePeriod_flag
            Aus_fun = @(idx) Aus_relsP.*ones(length(idx),1);
        end
        else
            % extract max dimension and check if matches number of oscillators
            [max_USwave_varlengths] = max(USwave_varlengths);
            idx_mUv_l = max_USwave_varlengths~=USwave_varlengths;
            if length(OSCindices)>max_USwave_varlengths
                error('size oscillators does not match dimensions USwaves')
            else
                % throw warning if less osccilators than dimensions
                if length(OSCindices)<max_USwave_varlengths
                    warning('USwaves more dimension than selected oscillators')
                end
                % check if all same dimensions if so do nothing
                if  any(USwave_varlengths-USwave_varlengths(1))
                    if any(USwave_varlengths(idx_mUv_l)-1)
                        error('length of US variables should either be the same for all or one should be bigger while others are one')
                    else
                        repmat_vals = max_USwave_varlengths-USwave_varlengths+1;
                    end
                    Dirus = repmat(Dirus,repmat_vals(1),1);
                    Aus = repmat(Aus,repmat_vals(2),1);
                    wus = repmat(wus,repmat_vals(3),1);
                    Phaseus = repmat(Phaseus,repmat_vals(4),1);
                    
                    
                end
                USwave_interm = @(t,idx) Dirus(idx,:).*Aus(idx,:).*sin(wus(idx).*t+Phaseus(idx,:));
                USwave_all = @(t,idx) USwave_interm(t,idx==OSCindices);
                if vibrParam_flag
                    Dirus_fun = @(idx) Dirus(idx==OSCindices,:);
                    Phaseus_fun = @(idx) Phaseus(idx==OSCindices,:);
                    Aus_fun = @(idx) Aus(intersect(idx,OSCindices),:);
                end
                if rel_singlePeriod_flag
                    Aus_fun = @(idx) Aus_relsP(intersect(idx,OSCindices),:);
                end
                
            end
        end
    end
else
    if AdaptVamp_flag
        SF = Aus_AVa/vmax_AVa;
    else
        SF = 1/wus;
    end
    USwave_all = @(t,idx) SF.*[vampx(idx).*sin(wus.*t+vphasex(idx)),vampy(idx).*sin(wus.*t+vphasey(idx)),vampz(idx).*sin(wus.*t+vphasez(idx))];
    if vibrParam_flag
        Dirus_fun = @(idx) SF.*[maxvx(idx),maxvy(idx),maxvz(idx)];
        Aus_fun = @(idx) SF.*norm([maxvx(idx),maxvy(idx),maxvz(idx)]);
        Phaseus_fun = @(idx) phasemaxdir(idx);
    end
    % velocity vectors are given als input v(r,t) = Re(v(r)exp(jwt))==>
    % v(r,t) =
    % [|vx(r)|cos(wt+thetax),|vy(r)|cos(wt+thetay),|vz(r)|cos(wt+thetaz)]
    % x(r,t) =  int(v(r,t),t) = 1/w*
    % [|vx(r)|sin(wt+thetax),|vy(r)|sin(wt+thetay),|vz(r)|sin(wt+thetaz)]

    
end



end