function [RMS,SNR] = getRMSnSNRptherm(VR,VRDOI,s,NoiseAmp_idx,Input,nos4l_flag,Param,varargin)
rng(s)
Settings = Input.Settings;
fus = Settings{find(strcmpi(Settings,'fus'))+1};
Window = 1000;
Tend = 0.025;
resUS = 20;
Alphafun = @(t,Tau,t0) double((t-t0)>=0).*((t-t0)./Tau).*exp(1-(t-t0)./Tau)+0;
Tau = 0.005;
AlphaDelay = 0;
Ifun = @(t) Alphafun(t,Tau,AlphaDelay);
meanI = 10;

if ~isempty(Param)
    Window = Param.fbandwidth;
    Tend = Param.Tend;
    resUS = Param.resUS;
    meanI = Param.meanI;
    Ifun = Param.Ifun;
    
end


Tend = ceil(fus*Tend)/fus; % recalculate Tend such that frequency spectrum contains fus
dt = (resUS*fus)^-1;
Tsim = 0:dt:Tend;
inputSignal = [Tsim;meanI.*Ifun(Tsim)];
mirror = zeros(size(VR,1),1);

if length(Tsim)~=size(VR,2)
    if length(Tsim)==size(VR,1)
        VR = VR';
    else
        error('wrong sizes')
    end
end
if length(Tsim)~=size(VRDOI,2)
    if length(Tsim)==size(VRDOI,1)
        VRDOI = VRDOI';
    else
        error('wrong sizes')
    end
end
wus = 2*pi*fus;
SF = 1/(wus);
NoiseAmps= [10^-2,10^-2*sqrt(2*Window),10^-2*sqrt(1/dt)];

% extract vibration direction to determine the mirror
if nos4l_flag
    DOIindice = 1;
    Dirus_way = varargin{1};
    CSource = varargin{2};
    CSink = varargin{3};
    Aus_way = varargin{4};
    Aus = varargin{5};
    k_Aus = varargin{6};
    Phaseus = varargin{7};
    Dirus = varargin{8};
    switch lower(Dirus_way)
        case 'dipole_dir'

            Dirus_DOI = Dirus(DOIindice,:);
            Dirus = CSource-CSink;
            Dirus = Dirus./vecnorm(Dirus);
            Dirus(DOIindice,:) = Dirus_DOI;
            
        case 'radial'
            
            Dirus_DOI = Dirus(DOIindice,:);
            Dirus = (CSource+CSink)/2;
            Dirus = Dirus./vecnorm(Dirus);
            Dirus(DOIindice,:) = Dirus_DOI;
        case 'input'
            if Display
                fprintf('\n Dirus not changed from input')
            end
        otherwise 
            error('false input Dirus_way')

    end
    switch lower(Aus_way)
        case 'kexp'
           
            if isrow(Aus)
                Aus = Aus';
            end
            Aus_DOI = Aus(DOIindice,:);
            kfun = @(x)  Aus_DOI*exp(-k_Aus*x);
            Rdiffdppos = vecnorm((CSource+CSink-CSource(DOIindice,:)-CSink(DOIindice,:))/2,2,2);
            Aus = kfun(Rdiffdppos);
            Aus(DOIindice,:) = Aus_DOI;            
            
        case 'kscale'

            if isrow(Aus)
                Aus = Aus';
            end
            Aus_DOI = Aus(DOIindice,:);
            Aus = Aus_DOI*k_Aus*size(CSource,1);
            Aus(DOIindice,:) = Aus_DOI;
            
        case 'normal'

        otherwise
            error('false input Aus_way')
    end
    USwave_varlengths = [size(Dirus,1), size(Aus,1), length(wus), size(Phaseus,1)];
    [max_USwave_varlengths, idx_mUv] = max(USwave_varlengths);
    idx_mUv_l = max_USwave_varlengths~=USwave_varlengths;
    if  any(USwave_varlengths-USwave_varlengths(1))
        if any(USwave_varlengths(idx_mUv_l)-1)
            error('length of US variables should either be the same for all or one should be bigger while others are one')
        else
            repmat_vals = max_USwave_varlengths-USwave_varlengths+1;
        end
        switch idx_mUv
            case 1
                Aus = repmat(Aus,repmat_vals(2),1);
                wus = repmat(wus,repmat_vals(3),1);
                Phaseus = repmat(Phaseus,repmat_vals(4),1);
            case 2
                Dirus = repmat(Dirus,repmat_vals(1),1);
                wus = repmat(wus,repmat_vals(3),1);
                Phaseus = repmat(Phaseus,repmat_vals(4),1);
            case 3
                Dirus = repmat(Dirus,repmat_vals(1),1);
                Aus = repmat(Aus,repmat_vals(2),1);
                Phaseus = repmat(Phaseus,repmat_vals(4),1);
            case 4
                Dirus = repmat(Dirus,repmat_vals(1),1);
                Aus = repmat(Aus,repmat_vals(2),1);
                wus = repmat(wus,repmat_vals(3),1);
        end
        
    end
    USwave_all = @(t,idx) Dirus(idx,:).*Aus(idx,:).*sin(wus(idx).*t+Phaseus(idx,:));
else
    vampx = varargin{1}; vampy = varargin{2}; vampz = varargin{3};
    vphasex = varargin{4}; vphasey = varargin{5}; vphasez = varargin{6};
    DOIindice = 1;
USwave_all = @(t,idx) SF.*[vampx(idx).*sin(wus.*t+vphasex(idx)),vampy(idx).*sin(wus.*t+vphasey(idx)),vampz(idx).*sin(wus.*t+vphasez(idx))];

end
% determine theta of USwave
tvals = linspace(0,1/fus,1e4)'; tvals = tvals(1:end-1);
USwave_sel = USwave_all(tvals,1);
USwave_sel = USwave_sel-mean(USwave_sel);
amp_vals = max(USwave_sel);
USwave_sel = USwave_sel(1:floor(length(tvals)/2),:);
USwave_t0 = USwave_sel(1,:);
dUSwave_t1t0= USwave_sel(2,:)-USwave_sel(1,:);

[~,cross_zero] = min(abs(USwave_sel));
theta_vals(1) = wus(DOIindice)*(1/(2*fus)*double(USwave_t0(1)>0)-tvals(cross_zero(1)));
theta_vals(2) = wus(DOIindice)*(1/(2*fus)*double(USwave_t0(2)>0)-tvals(cross_zero(2)));
theta_vals(3) = wus(DOIindice)*(1/(2*fus)*double(USwave_t0(3)>0)-tvals(cross_zero(3)));
theta_vals(cross_zero==1) = 0; theta_vals(cross_zero==1 & dUSwave_t1t0<0) = pi;
% USwave_t0 = USwave_sel(t0);USwave_t1 = USwave_sel(t1);
% theta_vals = atan2(sin(wus.*t1).*ones(1,3),(USwave_t1./USwave_t0-cos(wus.*t1)));
% % nanmean is for special case when theta is zero ==> sin(theta) is 0 gives
% % devision by zero = nan
% amp_vals = nanmean(vertcat(abs(USwave_t0./sin(theta_vals)),abs(USwave_t1./sin(wus.*t1+theta_vals))));

% circular mean of different directions
avg_vect = nansum(amp_vals.*exp(complex(0,theta_vals)))./nansum(amp_vals);
theta_global = atan2(imag(avg_vect),real(avg_vect));


VRnoise =  VR-VRDOI;
POIsidx = 1:size(VR,1);
noise = NoiseAmps(NoiseAmp_idx)*randn(size(VR));
VRptherm = VR+noise; VRptherm_noise = VRnoise+noise;

[reconsSignalVR,powerSignalVR,timeReconS] = calcFourierGetsignal(...
    VRptherm,Tsim,POIsidx,fus,Window,inputSignal,mirror,theta_global,0);
[reconsSignalVRDOI,powerSignalVRDOI] = calcFourierGetsignal(...
    VRDOI,Tsim,POIsidx,fus,Window,inputSignal,mirror,theta_global,0);
[reconsSignalVRnoise,powerSignalVRnoise] = calcFourierGetsignal(...
    VRptherm_noise,Tsim,POIsidx,fus,Window,inputSignal,mirror,theta_global,0);
for iPOI=POIsidx
    SNR(iPOI) = 10*log10(powerSignalVRDOI(iPOI)/powerSignalVRnoise(iPOI));
    [RMS(iPOI)] = ...
        calcRMS(reconsSignalVR{iPOI},inputSignal,timeReconS{iPOI},mirror(iPOI),0);
end
end