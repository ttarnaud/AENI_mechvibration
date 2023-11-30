function [VRMAT,VRMATosc,VRMATnoise]=SimDipoleOscPC(Tend,CSource,CSink,AmpI,USwave,USperiod,idxOscDip,varargin)

%default settings
Validation_flag = 0;
display_flag = 1;
sigma = 0.33; %S/m
RSphere = 0.07; %m
resUS = 100;
POIs = [];      %declaration POIs
PlotSS_flag = 0;
SolutionType = '3Sphere';
Include_frequency_flag = 1;
scale_flag = 1;
% other sphere thicknesses
tSkull = 0.005;
tScalp = 0.007;
tAir = 0.005;

% Declare current in dipoles
idxI = ones(size(CSource,1),1);
Istandard = 1;          % standard is 1 µA
I = Istandard*ones(size(CSource,1),1);

if iscell(varargin) && length(varargin) == 1 && iscell(varargin{1})
    varargin = varargin{1};
end

if mod(length(varargin)+1,2)
    if ~isempty(varargin)
        if any(strcmpi(varargin,'display'))
            display_flag = varargin{find(strcmpi(varargin,'display'))+1};
        end
        if any(strcmpi(varargin,'RSphere'))
            RSphere = varargin{find(strcmpi(varargin,'RSphere'))+1};
        end
        if any(strcmpi(varargin,'sigma'))
            sigma = varargin{find(strcmpi(varargin,'sigma'))+1};
        end
        if any(strcmpi(varargin,'resUS'))
            resUS = varargin{find(strcmpi(varargin,'resUS'))+1};
        end
        if any(strcmpi(varargin,'idxI'))
            idxI = varargin{find(strcmpi(varargin,'idxI'))+1};
        end
        if any(strcmpi(varargin,'POI'))
            POIs = varargin{find(strcmpi(varargin,'POI'))+1};
        end
        if any(strcmpi(varargin,'SolutionType'))
            SolutionType = varargin{find(strcmpi(varargin,'SolutionType'))+1};
        end
        if any(strcmpi(varargin,'fDependence'))
            Include_frequency_flag = varargin{find(strcmpi(varargin,'fDependence'))+1};
        end
        if any(strcmpi(varargin,'tSkull'))
            tSkull = varargin{find(strcmpi(varargin,'tSkull'))+1};
        end
        if any(strcmpi(varargin,'tScalp'))
            tScalp = varargin{find(strcmpi(varargin,'tScalp'))+1};
        end
        if any(strcmpi(varargin,'tAir'))
            tAir = varargin{find(strcmpi(varargin,'tAir'))+1};
        end
        if any(strcmpi(varargin,'scale'))
            scale_flag = varargin{find(strcmpi(varargin,'scale'))+1};
        end
        if any(strcmpi(varargin,'Validation'))
            Validation_flag = varargin{find(strcmpi(varargin,'Validation'))+1};
        end
        
    end
else
    error('incorrect input')
end


% declare oscilating dipoles
CSourceinit = CSource;
CSinkinit = CSink;


% simulate first single oscilation cycle
dt = USperiod/resUS;
TsimUS = 0:dt:USperiod-dt;
if isempty(TsimUS)
    TsimUS = 1;
    disp('resolution not high enough for US modeling')
end
Tsim = 0:dt:Tend;
f = ceil(length(Tsim)/length(TsimUS));

% generate I matrix
if isa(AmpI,'function_handle')
    I = Istandard*ones(size(CSource,1),length(Tsim));  % elongate making it a matrix
    I(logical(idxI),:) = AmpI(Tsim); % every point of idxI is now changed by evaluation of Ampi at time t e Tsim
elseif (size(AmpI,1)==sum(idxI) || length(AmpI)==1) && size(AmpI,2) == 1 % this is for amplifying certain dipoles with constant value for whole sim time
    I(logical(idxI)) = AmpI;
    I = repmat(I,1,length(Tsim));
elseif isequal(size(AmpI),[size(CSource,1),length(Tsim)])
    I = AmpI;
elseif  (size(AmpI,1)==sum(idxI) || size(AmpI,1)==1) && size(AmpI,2) == length(Tsim)
    I = Istandard*ones(size(CSource,1),length(Tsim));
    I(logical(idxI),:) = AmpI;
else
    error('false ApI input')
end



switch Include_frequency_flag
    case 0
        VRMAT = zeros(size(POIs,1),length(Tsim));
        VRMATosc = zeros(size(POIs,1),length(Tsim));
        VRMATnoise = zeros(size(POIs,1),length(Tsim));
        parfor idp = 1:size(CSource,1) % loop over dipoles
            if norm(CSource(idp,:))>=RSphere || norm(CSink(idp,:))>=RSphere
                error('dipole not within brain region')
            end
            VRVibration = zeros(size(POIs,1),length(TsimUS));
            Inputdipole = struct();
            if any(idp==idxOscDip) % check if the one selected is vibrating
                for itime=1:length(TsimUS) % calculation vibration effect
                    t = TsimUS(itime);
                    Inputdipole.CSource = CSourceinit(idp,:)+USwave(t);
                    Inputdipole.CSink = CSinkinit(idp,:)+USwave(t);
                    if ~isempty(POIs)
                        if size(POIs,2)==3
                            VRVibration(:,itime) = PotentialSingleSource(Inputdipole,abs(1*10^-6),sigma,POIs(:,1),POIs(:,2),POIs(:,3),RSphere,'',SolutionType,varargin);
                        else
                            error('size POIs incorrect')
                        end
                    end
                end
            else % if not oscilator calculated potential and elongate
                Inputdipole.CSource = CSourceinit(idp,:);
                Inputdipole.CSink = CSinkinit(idp,:);
                if ~isempty(POIs)
                    if size(POIs,2)==3
                        VRPOIs = PotentialSingleSource(Inputdipole,abs(1*10^-6),sigma,POIs(:,1),POIs(:,2),POIs(:,3),RSphere,'',SolutionType,varargin);
                        VRVibration = repmat(VRPOIs,1,length(TsimUS));
                    else
                        error('size POIs incorrect')
                    end
                end
            end
            VRVibration = repmat(VRVibration,1,f);
            VRMAT = VRMAT + I(idp,:).*VRVibration(:,1:length(I(idp,:)));
            if any(idp==idxOscDip)
                VRMATosc = VRMATosc + I(idp,:).*VRVibration(:,1:length(I(idp,:)));
            else
                VRMATnoise = VRMATnoise + I(idp,:).*VRVibration(:,1:length(I(idp,:)));
            end
        end
        VRMAT = VRMAT.*1e6; %outcome in µV
        VRMATosc = VRMATosc*1e6;
        VRMATnoise = VRMATnoise*1e6;
    case 1
        VRMAT = zeros(size(POIs,1),length(Tsim));
        VRMATosc = zeros(size(POIs,1),length(Tsim));
        VRMATnoise = zeros(size(POIs,1),length(Tsim));
        relTOL = 1e-10;
        absTOL = 1e-13;
        Fs = 1/(Tsim(2)-Tsim(1));  % extracting sampling frequency
        N = length(Tsim);           % extracting signal length
        RScalp = RSphere+tScalp+tSkull;  % distance to outerside of head boundary scalp/air
        RSkull = RSphere+tSkull; % distance to boundary skull/scalp
        RBrain = RSphere; % distance to boundary brain/skull
        nPOI = POIs./vecnorm(POIs,2,2); % normalized vector POIs [Nx3]
        Dataload = load('DataGabriels.mat');
        sigmaSkullfun = @(x) interp1(Dataload.Data2.Skull(:,1),Dataload.Data2.Skull(:,3),x);
        epsilonSkullfun = @(x) interp1(Dataload.Data2.Skull(:,1),Dataload.Data2.Skull(:,2),x);
        sigmaGMfun = @(x) interp1(Dataload.Data2.GreyMatter(:,1),Dataload.Data2.GreyMatter(:,3),x);
        epsilonGMfun = @(x) interp1(Dataload.Data2.GreyMatter(:,1),Dataload.Data2.GreyMatter(:,2),x);
        if strcmpi(SolutionType,'3Sphere')
            f1 = RBrain/RScalp; % normalized distance brain
            f2 = RSkull/RScalp; % normalized distance skull
        elseif strcmpi(SolutionType,'4Sphere')
            RAir = RSphere+tScalp+tSkull+tAir;
            eps0 = 8.85*10^-12;
            sigmaScalpfun = @(x) interp1(Dataload.Data2.Scalp(:,1),Dataload.Data2.Scalp(:,3),x);
            epsilonScalpfun = @(x) interp1(Dataload.Data2.Scalp(:,1),Dataload.Data2.Scalp(:,2),x);
        end
        parfor idp = 1:size(CSource,1) % loop over dipoles
            if norm(CSource(idp,:))>=RSphere || norm(CSink(idp,:))>=RSphere
                error('dipole not within brain region')
            end
            iter_flag = true; % flag for while loop
            iSum = 0;
            iSum_lim = 1e5;
            PosCDInit = (CSourceinit(idp,:)+CSinkinit(idp,:))./2;
            if any(PosCDInit)
                b0 = norm(PosCDInit);
            else
                b0 = norm(CSourceinit(idp,:));
            end
            if display_flag
                fprintf('Dipole %d:  \n',idp)
            end
            while (iter_flag || iSum<10) && iSum<iSum_lim
                iSum = iSum+1;
                fitime = zeros(size(POIs,1),length(TsimUS));
                if any(idp==idxOscDip) % check if the one selected is vibrating
                    Pflag = 0;
                    for itime=1:length(TsimUS) % calculation vibration effect
                        t = TsimUS(itime);
                        CSource_Sel = CSourceinit(idp,:) + USwave(t);
                        CSink_Sel = CSinkinit(idp,:) + USwave(t);
                        fitime(:,itime) = calcfitime(CSource_Sel,CSink_Sel,iSum,nPOI,b0);
                        if display_flag && iSum == 1
                            if round(itime/length(TsimUS)*100,0) >= Pflag
                                fprintf('\t OScillation: %d%%\n',Pflag)
                                Pflag = Pflag+10;
                            end
                        end
                    end
                else
                    CSource_Sel = CSourceinit(idp,:);
                    CSink_Sel = CSinkinit(idp,:);
                    fitime = calcfitime(CSource_Sel,CSink_Sel,iSum,nPOI,b0);
                    fitime = repmat(fitime,1,length(TsimUS));
                end
                fitime = repmat(fitime,1,f);
                fitime = I(idp,:).*fitime(:,1:length(I(idp,:))); % elongate and multiply with correct amplitude
                Fomega = fft(fitime,N,2);
                freqs = 0:Fs/N:Fs/2; f_neg=-Fs/N:-Fs/N:-Fs/2;% calculating negative components
                freqs = [freqs,fliplr(f_neg)];  % make one array
                if strcmpi(SolutionType,'3sphere')
                    Gomega = calcGomega(freqs,iSum,RScalp,f1,f2,sigmaSkullfun,sigmaGMfun,epsilonSkullfun,epsilonGMfun,b0);
                elseif strcmpi(SolutionType,'4Sphere')
                    sigmaBrain = complex(sigmaGMfun(abs(freqs)),freqs.*epsilonGMfun(abs(freqs)).*eps0);
                    sigmaSkull = complex(sigmaSkullfun(abs(freqs)),freqs.*epsilonSkullfun(abs(freqs)).*eps0);
                    sigmaScalp = complex(sigmaScalpfun(abs(freqs)),freqs.*epsilonScalpfun(abs(freqs)).*eps0);
                    sigmaAir = complex(0.*ones(size(freqs)),freqs.*ones(size(freqs)).*eps0);
                    if Validation_flag
                        sigmaBrain = sigma(1).*ones(size(freqs));
                        sigmaSkull = sigma(2).*ones(size(freqs));
                        sigmaScalp = sigma(3).*ones(size(freqs));
                        sigmaAir = sigma(4).*ones(size(freqs));
                    end
                    if scale_flag
                        Gomega = GetGi(RAir*ones(size(POIs,1),1),RBrain,RSkull,RScalp,RAir,sigmaBrain,sigmaSkull,sigmaScalp,sigmaAir,iSum,b0);
                    else
                        Gomega = GetGi(vecnorm(POIs,2,2),RBrain,RSkull,RScalp,RAir,sigmaBrain,sigmaSkull,sigmaScalp,sigmaAir,iSum,b0);
                    end
                else
                    error('no frequency dependence combined with this type of solution')
                end
                VRiSumf = Gomega.*Fomega;
                VRiSumt = ifft(VRiSumf,N,2,'symmetric');
                if any(isnan(VRiSumt(:)))
                    disp(['nan at ',num2str(iSum)])
                end
                if iSum ~= 1
                    VRMATssOLD = VRssMAT;
                    VRssMAT = VRssMAT+VRiSumt;
                    diffVRssMAT = abs(VRssMAT-VRMATssOLD);
                    if (max(diffVRssMAT(:)./VRMATssOLD(:)) < relTOL && max(diffVRssMAT(:)) <absTOL) || max(diffVRssMAT(:)) == 0
                        iter_flag = false;
                    end
                    
                    if iSum == iSum_lim
                        disp('max iterations reached')
                        disp(['max add value',num2str(max(abs(VRssMAT(:)-VRMATssOLD(:))./VRMATssOLD(:)))])
                        break
                    end
                    if ~logical(mod(iSum,10)) && display_flag
                        disp(['Legendre P needed: ',num2str(iSum)])
                    end
                else
                    VRssMAT = VRiSumt;
                end
            end
            VRMAT = VRMAT+VRssMAT*1e6; % output is in µV
            if any(idp==idxOscDip)
                VRMATosc = VRMATosc+VRssMAT*1e6; % output is in µV
            else
                VRMATnoise = VRMATnoise+VRssMAT*1e6; % output is in µV
            end
        end
    otherwise
        error('wrong frequency dependence input')
end

end
