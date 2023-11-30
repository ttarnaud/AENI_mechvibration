function [VRMAT,VRMATDOI,VRMATosc,VRMATnoise] = SimDipoleOsc_DB_OLD(Tend,CSource,CSink,Iarray,USwave,USperiod,idxOscDip,varargin)

vib_interp_type = 'linear';
Display = 1;
%% change default variables
if iscell(varargin) && length(varargin) == 1 && iscell(varargin{1})
    varargin = varargin{1};
end
if mod(length(varargin)+1,2)
    if ~isempty(varargin)
        if any(strcmpi(varargin,'vib_interp_type'))
            vib_interp_type = varargin{find(strcmpi(varargin,'vib_interp_type'))+1};
        end
        if any(strcmpi(varargin,'Display'))
            Display = varargin{find(strcmpi(varargin,'Display'))+1};
        end
    end
end

fus = 1/USperiod;
wus = 2*pi*fus;
dp_pos = (CSource+CSink)/2;
        



try
    SolutionType = varargin{find(strcmpi(varargin,'SolutionType'))+1};
    sigma = varargin{find(strcmpi(varargin,'sigma'))+1};
    idxDOIindice = varargin{find(strcmpi(varargin,'DOI'))+1};
    Include_frequency_flag = varargin{find(strcmpi(varargin,'fDependence'))+1};
    POIs = varargin{find(strcmpi(varargin,'POI'))+1};
    RSphere = varargin{find(strcmpi(varargin,'RSphere'))+1};
    resUS = varargin{find(strcmpi(varargin,'resUS'))+1};
    tSkull = varargin{find(strcmpi(varargin,'tSkull'))+1};
    tScalp = varargin{find(strcmpi(varargin,'tScalp'))+1};
    tAir = varargin{find(strcmpi(varargin,'tAir'))+1};
catch
    error('at least provide DOI indice,\n fDependence, POI, RSphere, resUS,\n tSkull, tScalp and tAir')
end


dt = (1/fus)/resUS;
TsimUS = 0:dt:(1/fus)-dt;
if isempty(TsimUS)
    TsimUS = 1;
    disp('resolution not high enough for US modeling')
end
Tsim = 0:dt:Tend;
f = ceil(length(Tsim)/length(TsimUS));

RPOIs = vecnorm(POIs,2,2); 
if any((RPOIs-RPOIs(1))>1e-6)
    error('POIs at different radii => no database available yet')
else
    RPOI = RPOIs(1);
end

switch RPOI
    case 0.082
    if (RPOI-(RSphere+tSkull+tScalp+tAir))>1e-6 || (fus-1e6)>1e-3 || Include_frequency_flag
        error('mismatch settings and database')
    end
    DB_prefix = 'Databasehuman_scalp_v2_1.00MHz_pos';
    pos_vals = 0:0.5:65;
    
    Aus_DB = 1e-6;%vibration amplitude was 1�m
    if ~contains(DB_prefix,'_v2')
        cf_DB = 1e-6; %database is recorded in �V => scale back to V
        cf_DB = cf_DB/10; % database recorded with Iarray 10�A
    else
        cf_DB = 1;
    end
    varargin{find(strcmpi(varargin,'Scale'))+1} = 0;
    otherwise
        error('no database yet')
end
% figure out position of dipoles
% sort in order to minimize database reloads
Rdps = round(vecnorm(dp_pos,2,2),6);
[Rdps_s,idx_s] = sort(Rdps);
Rdps_u = unique(Rdps_s);
% sort other inputs
Iarray = Iarray(idx_s,:);
USwave_s =@(t,i) USwave(t,idx_s(i));
dp_pos = dp_pos(idx_s,:);
% new positions of oscDipoles
[~,~,idxOscDip] = intersect(idxOscDip,idx_s,'stable');
CSource = CSource(idx_s,:);
CSink = CSink(idx_s,:);

if strcmpi(vib_interp_type,'linear')
    try
        Dirus = varargin{find(strcmpi(varargin,'Dirus'))+1};
        Phaseus = varargin{find(strcmpi(varargin,'Phaseus'))+1};
    catch
        error('if linear mode include Dirus and Phaseus')
    end
    Dirus_s =@(idx) Dirus(idx_s(idx))/norm(Dirus(idx_s(idx)));
    Aus_s = @(idx) norm(Dirus(idx_s(idx)))/Aus_DB; % Database recorded at vibration amplitude = |v|/wus =  1e-6
    Phaseus_s =@(idx) Phaseus(idx_s(idx));
end

VRMAT = zeros(size(POIs,1),length(Tsim));
VRMATDOI = zeros(size(POIs,1),length(Tsim));
VRMATosc = zeros(size(POIs,1),length(Tsim));
VRMATnoise = zeros(size(POIs,1),length(Tsim));
for iru = 1:length(Rdps_u)
    Rdp = Rdps_u(iru);
    % find closest lower and closest higher
    pos_vals_intm = pos_vals-Rdp;
    pos_vals_lower = pos_vals_intm(pos_vals_intm<=0);
    pos_vals_higher = pos_vals_intm(pos_vals_intm>=0);
    Rdp_cl = max(pos_vals_lower)+Rdp;
    Rdp_ch = min(pos_vals_higher)+Rdp;
    % check if we need to (re)load a database
    if iru~=1
        load_flag = ~(Rdp_cl==Rdp_cl_loaded && Rdp_ch==Rdp_ch_loaded);
    else
        load_flag = 1;
    end
    if load_flag
        if Rdp_cl==Rdp_ch
            tDB_flag = 0; %two database flag
            DB_sufix = [sprintf('%0.2f',Rdp_cl),'mm.mat'];
            Database = load([DB_prefix,DB_sufix]);
        else
            tDB_flag = 1;
            DB_sufix_cl = [sprintf('%0.2f',Rdp_cl),'mm.mat'];
            DB_sufix_ch = [sprintf('%0.2f',Rdp_ch),'mm.mat'];
            Database = load([DB_prefix,DB_sufix_cl]);
            Database_ch = load([DB_prefix,DB_sufix_ch]);
        end
        Rdp_cl_loaded = Rdp_cl;
        Rdp_ch_loaded = Rdp_ch;
    end
    theta_POIs = acos(Database.POIs(:,3)/RPOI);
    
    dt_DB = fus^(-1)./(size(Database.DataBase,3)-1);
    TsimUS_DB = 0:dt_DB:fus^(-1);
    idx_dpOI = find(Rdps_s==Rdp);
    
    for idp=1:length(idx_dpOI)
        VRMAT_dpOI_vibr = [];
        Inputdipole = struct();
        dpOI = dp_pos(idx_dpOI(idp),:);
        IarrayOI = Iarray(idx_dpOI(idp),:);
        USwaveOI = USwave_s(TsimUS',idx_dpOI(idp));
        if any(idx_dpOI(idp)==idxOscDip)
            switch lower(vib_interp_type)
                case 'linear'
                    % linear interpolation based on amplitude and phase of
                    % field oscillation
                    VR_DB_dpOI_POIs = zeros(size(POIs,1),1);
                    phi_dpOI = zeros(size(POIs,1),1);
                    theta_dpOI = zeros(size(POIs,1),1);
                    % linear interpolation (dominant direction is given in
                    % input together with phase
                    DirusOI = Dirus_s(idx_dpOI(idp));
                    PhaseusOI = Phaseus_s(idx_dpOI(idp));
                    AusOI = Aus_s(idx_dpOI(idp));
                    
                    % determine angle between dpOI and POIs
                    if any(dpOI)
                        thetaOI = acos(1./(norm(dpOI)*RPOI).*dpOI*POIs');
                    else
                        thetaOI = acos(1./(norm(CSource(idx_dpOI(idp),:))*RPOI).*CSource(idx_dpOI(idp),:)*POIs');
                    end
                    thetaOI = round(thetaOI,4);
                    theta_intm = bsxfun(@minus,theta_POIs,thetaOI);
                    theta_cl = theta_intm; theta_cl(theta_cl>0) = -inf;
                    [theta_cl,idx_cl] = max(theta_cl);
                    theta_cl = theta_cl+thetaOI;
                    theta_ch = theta_intm; theta_ch(theta_ch<0) = inf;
                    [theta_ch,idx_ch] = min(theta_ch);
                    theta_ch = theta_ch+thetaOI;
                    
                    % interpolate to correct vib direction
                    % determine rotation matrix in order to be the dipoleOI on
                    % z-axis and POI in YZ plane
                    ndpDB_POIDB = cross(Database.POIs(end-1,:),Database.selpos,2);
                    ndpDB_POIDB = ndpDB_POIDB/norm(ndpDB_POIDB);
                    if any(isnan(ndpDB_POIDB)) && Display
                        warning('normal could not be determined, +x-axis is taken')
                        ndpDB_POIDB = [1,0,0];
                    end
                    % loop over different POIs
                    for iPOI = 1:size(POIs,1)
                        if ~any(dpOI)
                            testOI = CSource(idx_dpOI(idp),:);
                        else
                            testOI = dpOI;
                        end
                        testOI = testOI/norm(testOI);
                        ndpOI_POI = cross(POIs(iPOI,:),testOI);
                        % if POI and dp on same vector
                        if ~any(ndpOI_POI)
                            %ndpOI_POI = cross(POIs(iPOI,:),CSource(idx_dpOI(idp),:));
                            %if ~any(ndpOI_POI)
                                ndpOI_POI = [0,0,1;0,1,0;-1,0,0]*testOI';
                                ndpOI_POI = ndpOI_POI';
%                                 if ~any(ndpOI_POI)
%                                     ndpOI_POI = [0,0,1;0,1,0;-1,0,0]*CSource(idx_dpOI(idp),:)';
%                                     ndpOI_POI = ndpOI_POI';
%                                 end
%                             end
                        end
                        ndpOI_POI = ndpOI_POI./norm(ndpOI_POI);
                        
                        % determine rotation matrix
                        %rotate first so plane determ by dp and POI is YZ
                        A_plane = calc_RM(ndpOI_POI,ndpDB_POIDB);
                        if any(isnan(A_plane(:)))
                            A_plane = eye(3);
                            if any(ndpOI_POI-ndpDB_POIDB)
                                A_plane = A_plane.*(-2*double((ndpOI_POI-ndpDB_POIDB)~=0)+1);
                            end
                            
                        end
                        dpOI_intm = (A_plane*testOI')';
                        dpOI_intm = dpOI_intm/norm(dpOI_intm);
                        % rotate dp is on z-axis (rotate around x)
                        A_x = calc_RM(dpOI_intm,[0,0,1]);
                        if any(isnan(A_x(:)))
                            A_x = eye(3);
                            if any(dpOI_intm-[0,0,1])
                                A_x = A_x.*(-2*double((dpOI_intm-[0,0,1])~=0)+1);
                            end                            
                        end
                        
                        A = A_x*A_plane;
                        if any(abs((A*testOI')'-[0,0,1])>1e-6)
                            error('rotation not good')
                        end
                        
                        
                        % rotate vibration direction
                        DirusOI_rot = (A*DirusOI')';
                        % determine phi and thetas
                        phi_dpOI(iPOI) = atan2(DirusOI_rot(2),DirusOI_rot(1));
                        theta_dpOI(iPOI) = acos(DirusOI_rot(3));
                    end
                        % interpolate closest lower POI closest lower dp
                        AMP_DB_clPOI_cldp = slerp2_sph(Database.phi_vdir,...
                            Database.theta_vdir,Database.amp(:,idx_cl),phi_dpOI,theta_dpOI,'plot',0);
                        AMP_DB_clPOI_cldp = diag(AMP_DB_clPOI_cldp)';
                        Phase_DB_clPOI_cldp = slerp2_sph(Database.phi_vdir,...
                            Database.theta_vdir,Database.phase(:,idx_cl),phi_dpOI,theta_dpOI,'plot',0);
                        Phase_DB_clPOI_cldp = diag(Phase_DB_clPOI_cldp)';
                        Stat_DB_clPOI_cldp = slerp2_sph(Database.phi_vdir,...
                            Database.theta_vdir,Database.statval(:,idx_cl),phi_dpOI,theta_dpOI,'plot',0);
                        Stat_DB_clPOI_cldp = diag(Stat_DB_clPOI_cldp)';
                        
                        if any(idx_cl~=idx_ch)
                        % interpolate closest higher POI closest lower dp
                        AMP_DB_chPOI_cldp = slerp2_sph(Database.phi_vdir,...
                            Database.theta_vdir,Database.amp(:,idx_ch),phi_dpOI,theta_dpOI,'plot',0);
                        AMP_DB_chPOI_cldp = diag(AMP_DB_chPOI_cldp)';
                        Phase_DB_chPOI_cldp = slerp2_sph(Database.phi_vdir,...
                            Database.theta_vdir,Database.phase(:,idx_ch),phi_dpOI,theta_dpOI,'plot',0);
                        Phase_DB_chPOI_cldp = diag(Phase_DB_chPOI_cldp)';
                        Stat_DB_chPOI_cldp = slerp2_sph(Database.phi_vdir,...
                            Database.theta_vdir,Database.statval(:,idx_ch),phi_dpOI,theta_dpOI,'plot',0);
                        Stat_DB_chPOI_cldp = diag(Stat_DB_chPOI_cldp)';
                        % interpolate to correct POI closest lower dp
                        AMP_DB_cldp = (AMP_DB_chPOI_cldp-AMP_DB_clPOI_cldp)./(theta_ch-theta_cl).*(thetaOI-theta_cl)+AMP_DB_clPOI_cldp;
                        Phase_DB_cldp = (Phase_DB_chPOI_cldp-Phase_DB_clPOI_cldp)./(theta_ch-theta_cl).*(thetaOI-theta_cl)+Phase_DB_clPOI_cldp;
                        Stat_DB_cldp = (Stat_DB_chPOI_cldp-Stat_DB_clPOI_cldp)./(theta_ch-theta_cl).*(thetaOI-theta_cl)+Stat_DB_clPOI_cldp;
                        end
                        AMP_DB_cldp(idx_cl==idx_ch) = AMP_DB_clPOI_cldp(idx_cl==idx_ch);
                        Phase_DB_cldp(idx_cl==idx_ch) = Phase_DB_clPOI_cldp(idx_cl==idx_ch);
                        Stat_DB_cldp(idx_cl==idx_ch) = Stat_DB_clPOI_cldp(idx_cl==idx_ch);
            
                        
                        
                        if tDB_flag
                            % interpolate closest lower POI closest higher dp
                            AMP_DB_clPOI_chdp = slerp2_sph(Database_ch.phi_vdir,...
                                Database_ch.theta_vdir,Database_ch.amp(:,idx_cl),phi_dpOI,theta_dpOI);
                            AMP_DB_clPOI_chdp = diag(AMP_DB_clPOI_chdp)';
                            Phase_DB_clPOI_chdp = slerp2_sph(Database_ch.phi_vdir,...
                                Database_ch.theta_vdir,Database_ch.phase(:,idx_cl),phi_dpOI,theta_dpOI);
                            Phase_DB_clPOI_chdp = diag(Phase_DB_clPOI_chdp)';
                            Stat_DB_clPOI_chdp = slerp2_sph(Database_ch.phi_vdir,...
                                Database_ch.theta_vdir,Database_ch.statval(:,idx_cl),phi_dpOI,theta_dpOI);
                            Stat_DB_clPOI_chdp = diag(Stat_DB_clPOI_chdp)';
                        
                            if any(idx_cl~=idx_ch)
                            % interpolate closest higher POI closest higher dp
                            AMP_DB_chPOI_chdp = slerp2_sph(Database_ch.phi_vdir,...
                                Database_ch.theta_vdir,Database_ch.amp(:,idx_ch),phi_dpOI,theta_dpOI);
                            AMP_DB_chPOI_chdp = diag(AMP_DB_chPOI_chdp)';
                            Phase_DB_chPOI_chdp = slerp2_sph(Database_ch.phi_vdir,...
                                Database_ch.theta_vdir,Database_ch.phase(:,idx_ch),phi_dpOI,theta_dpOI);
                            Phase_DB_chPOI_chdp = diag(Phase_DB_chPOI_chdp)';
                            Stat_DB_chPOI_chdp = slerp2_sph(Database_ch.phi_vdir,...
                                Database_ch.theta_vdir,Database_ch.statval(:,idx_ch),phi_dpOI,theta_dpOI);
                            Stat_DB_chPOI_chdp = diag(Stat_DB_chPOI_chdp)';
                            % interpolate to correct POI closest higher dp
                            AMP_DB_chdp = (AMP_DB_chPOI_chdp-AMP_DB_clPOI_chdp)./(theta_ch-theta_cl).*(thetaOI-theta_cl)+AMP_DB_clPOI_chdp;
                            Phase_DB_chdp = (Phase_DB_chPOI_chdp-Phase_DB_clPOI_chdp)./(theta_ch-theta_cl).*(thetaOI-theta_cl)+Phase_DB_clPOI_chdp;
                            Stat_DB_chdp = (Stat_DB_chPOI_chdp-Stat_DB_clPOI_chdp)./(theta_ch-theta_cl).*(thetaOI-theta_cl)+Stat_DB_clPOI_chdp;
                            end
                            AMP_DB_chdp(idx_cl==idx_ch) = AMP_DB_clPOI_chdp(idx_cl==idx_ch);
                            Phase_DB_chdp(idx_cl==idx_ch) = Phase_DB_clPOI_chdp(idx_cl==idx_ch);
                            Stat_DB_chdp(idx_cl==idx_ch) = Stat_DB_clPOI_chdp(idx_cl==idx_ch);
                        
                            % interpolate to correct dp
                            AMP_DB = (AMP_DB_chdp-AMP_DB_cldp)./(Rdp_ch-Rdp_cl).*Rdp+AMP_DB_cldp;
                            Phase_DB = (Phase_DB_chdp-Phase_DB_cldp)./(Rdp_ch-Rdp_cl).*Rdp+Phase_DB_cldp;
                            Stat_DB = (Stat_DB_chdp-Stat_DB_cldp)./(Rdp_ch-Rdp_cl).*Rdp+Stat_DB_cldp;
                        else
                            AMP_DB = AMP_DB_cldp;
                            Phase_DB = Phase_DB_cldp;
                            Stat_DB = Stat_DB_cldp;
                        end
                        % interpolate for difference in time points
                        VRMAT_dpOI_vibr = AusOI*AMP_DB'.*sin(wus*TsimUS-Phase_DB'+PhaseusOI)+Stat_DB'; 

                case 'linear_tp'
                    % linear interpolation for each time point
                    VR_DB_dpOI_POIs = zeros(size(POIs,1),length(TsimUS));
                    % linear interpolation (dominant direction is given in
                    % input together with phase
                    DirusOI = Dirus_s(idx_dpOI(idp));
                    PhaseusOI = Phaseus_s(idx_dpOI(idp));
                    AusOI = Aus_s(idx_dpOI(idp));
                    
                    % determine angle between dpOI and POIs
                    if any(dpOI)
                        thetaOI = acos(1./(norm(dpOI)*RPOI).*dpOI*POIs');
                    else
                        thetaOI = acos(1./(norm(CSource(idx_dpOI(idp),:))*RPOI).*CSource(idx_dpOI(idp),:)*POIs');
                    end
                    thetaOI = round(thetaOI,4);
                    theta_intm = bsxfun(@minus,theta_POIs,thetaOI);
                    theta_cl = theta_intm; theta_cl(theta_cl>0) = -inf;
                    [theta_cl,idx_cl] = max(theta_cl);
                    theta_cl = theta_cl+thetaOI;
                    theta_ch = theta_intm; theta_ch(theta_ch<0) = inf;
                    [theta_ch,idx_ch] = min(theta_ch);
                    theta_ch = theta_ch+thetaOI;
                    
                    % interpolate to correct vib direction
                    % determine rotation matrix in order to be the dipoleOI on
                    % z-axis and POI in YZ plane
                    ndpDB_POIDB = cross(Database.POIs(end-1,:),Database.selpos,2);
                    ndpDB_POIDB = ndpDB_POIDB/norm(ndpDB_POIDB);
                    if any(isnan(ndpDB_POIDB))
                        warning('normal could not be determined, +x-axis is taken')
                        ndpDB_POIDB = [1,0,0];
                    end
                    % loop over different POIs
                    for iPOI = 1:size(POIs,1)
                        if ~any(dpOI)
                            testOI = CSource(idx_dpOI(idp),:);
                        else
                            testOI = dpOI;
                        end
                        testOI = testOI/norm(testOI);
                        ndpOI_POI = cross(POIs(iPOI,:),testOI);
                        % if POI and dp on same vector
                        if ~any(ndpOI_POI)
                            %ndpOI_POI = cross(POIs(iPOI,:),CSource(idx_dpOI(idp),:));
                            %if ~any(ndpOI_POI)
                                ndpOI_POI = [0,0,1;0,1,0;-1,0,0]*testOI';
                                ndpOI_POI = ndpOI_POI';
%                                 if ~any(ndpOI_POI)
%                                     ndpOI_POI = [0,0,1;0,1,0;-1,0,0]*CSource(idx_dpOI(idp),:)';
%                                     ndpOI_POI = ndpOI_POI';
%                                 end
%                             end
                        end
                        ndpOI_POI = ndpOI_POI./norm(ndpOI_POI);
                        
                        % determine rotation matrix
                        %rotate first so plane determ by dp and POI is YZ
                        A_plane = calc_RM(ndpOI_POI,ndpDB_POIDB);
                        if any(isnan(A_plane(:)))
                            A_plane = eye(3);
                            if any(ndpOI_POI-ndpDB_POIDB)
                                A_plane = A_plane.*(-2*double((ndpOI_POI-ndpDB_POIDB)~=0)+1);
                            end
                            
                        end
                        dpOI_intm = (A_plane*testOI')';
                        dpOI_intm = dpOI_intm/norm(dpOI_intm);
                        % rotate dp is on z-axis (rotate around x)
                        A_x = calc_RM(dpOI_intm,[0,0,1]);
                        if any(isnan(A_x(:)))
                            A_x = eye(3);
                            if any(dpOI_intm-[0,0,1])
                                A_x = A_x.*(-2*double((dpOI_intm-[0,0,1])~=0)+1);
                            end                            
                        end
                        
                        A = A_x*A_plane;
                        if any(abs((A*testOI')'-[0,0,1])>1e-6)
                            error('rotation not good')
                        end
                        
                        
                        % rotate vibration direction
                        DirusOI_rot = (A*DirusOI')';
                        % determine phi and thetas
                        phi_dpOI = atan2(DirusOI_rot(2),DirusOI_rot(1));
                        theta_dpOI = acos(DirusOI_rot(3));
                        % interpolate closest lower POI closest lower dp
                        VR_DB_clPOI_cldp = slerp2_sph(Database.phi_vdir,...
                            Database.theta_vdir,squeeze(Database.DataBase(:,idx_cl(iPOI),:)),phi_dpOI,theta_dpOI,'plot',0);
                        if idx_cl(iPOI)~= idx_ch(iPOI)
                        % interpolate closest higher POI closest lower dp
                        VR_DB_chPOI_cldp = slerp2_sph(Database.phi_vdir,...
                            Database.theta_vdir,squeeze(Database.DataBase(:,idx_ch(iPOI),:)),phi_dpOI,theta_dpOI,'plot',0);
                        % interpolate to correct POI closest lower dp
                        VR_DB_cldp = (VR_DB_chPOI_cldp-VR_DB_clPOI_cldp)./(theta_ch(iPOI)-theta_cl(iPOI)).*(thetaOI(iPOI)-theta_cl(iPOI))+VR_DB_clPOI_cldp;
                        else
                            VR_DB_cldp = VR_DB_clPOI_cldp;
                        end
                        if tDB_flag
                            % interpolate closest lower POI closest higher dp
                            VR_DB_clPOI_chdp = slerp2_sph(Database_ch.phi_vdir,...
                                Database_ch.theta_vdir,squeeze(Database_ch.DataBase(:,idx_cl(iPOI),:)),phi_dpOI,theta_dpOI);
                            if idx_cl(iPOI)~= idx_ch(iPOI)
                            % interpolate closest higher POI closest higher dp
                            VR_DB_chPOI_chdp = slerp2_sph(Database_ch.phi_vdir,...
                                Database_ch.theta_vdir,squeeze(Database_ch.DataBase(:,idx_ch(iPOI),:)),phi_dpOI,theta_dpOI);
                            % interpolate to correct POI closest higher dp
                            VR_DB_chdp = (VR_DB_chPOI_chdp-VR_DB_clPOI_chdp)./(theta_ch(iPOI)-theta_cl(iPOI)).*(thetaOI(iPOI)-theta_cl(iPOI))+VR_DB_clPOI_chdp;
                            else
                                VR_DB_chdp = VR_DB_clPOI_chdp;
                            end
                            % interpolate to correct dp
                            VR_DB = (VR_DB_chdp-VR_DB_cldp)./(Rdp_ch-Rdp_cl).*Rdp+VR_DB_cldp;
                        else
                            VR_DB = VR_DB_cldp;
                        end
                        % interpolate for difference in time points
                        VR_DB_dpOI_POIs(iPOI,:) = interp1(TsimUS_DB,VR_DB,TsimUS);
                    end
                    % include phase difference
                    nri_shift = floor(PhaseusOI/(2*pi*fus)/dt);
                    VRMAT_dpOI_vibr = circshift(VR_DB_dpOI_POIs,-nri_shift,2);
                    VRMAT_dpOI_vibr = cf_DB*VRMAT_dpOI_vibr;
                    VRMAT_dpOI_vibr = AusOI*(VRMAT_dpOI_vibr-mean(VRMAT_dpOI_vibr))+mean(VRMAT_dpOI_vibr);
                    
                case {'complete'}
                    error('not constructed yet')
                otherwise
                    error('type not defined')
            end

            
        else
            Inputdipole.CSource = CSource(idx_dpOI(idp),:);
            Inputdipole.CSink = CSink(idx_dpOI(idp),:);
            if ~isempty(POIs)
                if size(POIs,2)==3
                    VRPOIs = PotentialSingleSource(Inputdipole,abs(1*10^-6),sigma,POIs(:,1),POIs(:,2),POIs(:,3),RSphere,'',SolutionType,varargin);
                    VRMAT_dpOI_vibr = repmat(VRPOIs,1,length(TsimUS));
                else
                    error('size POIs incorrect')
                end
            end
        end

        
        VRMAT_dpOI_vibr = repmat(VRMAT_dpOI_vibr,1,f);
        VRMAT = VRMAT + IarrayOI.*VRMAT_dpOI_vibr(:,1:length(IarrayOI));
        if any(idx_s(idx_dpOI(idp)) == idxDOIindice)
            VRMATDOI = VRMATDOI + IarrayOI.*VRMAT_dpOI_vibr(:,1:length(IarrayOI));
        end
        if any(idx_dpOI(idp)==idxOscDip)
            VRMATosc = VRMATosc + IarrayOI.*VRMAT_dpOI_vibr(:,1:length(IarrayOI));
        else
            VRMATnoise = VRMATnoise + IarrayOI.*VRMAT_dpOI_vibr(:,1:length(IarrayOI));
        end
    end
end
VRMAT = VRMAT.*1e6; % outcome is given in �V
VRMATDOI = VRMATDOI.*1e6;
VRMATosc = VRMATosc*1e6;
VRMATnoise = VRMATnoise*1e6;

end