% Create results figure 2 for paper
% results from ampus study (change in spatial constant of ultrasonic field)
close all; clear all; clc
folder_rawresult = 'D:\no backup\EEGUS\HPC_files\1231\Crop';
filename_rawresult = 'ResultsCrop_cEG_Settingset1_Amp_fkApdpset_NCfixed_nos4l_v2_12-30-20_1825';
data_rr =  load(fullfile(folder_rawresult,filename_rawresult));


Alphafun = @(t,Tau,t0) double((t-t0)>=0).*((t-t0)./Tau).*exp(1-(t-t0)./Tau)+0; %alpha function
Tau = 0.005;                     % [s]
AlphaDelay = 0;                 %delay in alpha function [s]
Ifun = @(t) Alphafun(t,Tau,AlphaDelay); %default function applied as current on DOI

[M,N] = size(data_rr.Outall);
Amps = nan(M,N);
kaus = nan(M,N);
for i =1:M
    for j =1:N
        kaus(i,j) = data_rr.Outall(i,j).Out.Param.k_Aus;
        Amps(i,j) = data_rr.Outall(i,j).Out.Param.Aus;
    end
end
kaus = 1./kaus(:,1).*1000;
Amps = unique(Amps)*1e6;
fprintf(['spatial constants k: ', repmat('%5.2f ', 1, length(kaus)),'\n'],kaus')
fprintf(['Amps: ', repmat('%5.2e ', 1, length(Amps)),'\n'],Amps')



%%
figure(10)
idp = 5;
plotmethod ='norm';
myvals = 0;
ikaus = 6;
styles = {'-','--'};

Tend = data_rr.Outall(ikaus,idp).Out.Param.Tend;
resUS = data_rr.Outall(ikaus,idp).Out.Param.resUS;
fus = data_rr.Outall(ikaus,idp).Out.Param.fus;
Tend = ceil(fus*Tend)/fus; % recalculate Tend such that frequency spectrum contains fus
dt = (resUS*fus)^-1;
t = 0:dt:Tend;

IOI = 10.*Ifun(t);
inputSignal = [t;IOI];

% plot alpha 
Tsim = 0:0.1:data_rr.trS_mat(ikaus,idp,end,end)*1000;
plot(Tsim,Ifun(Tsim/1000),'color',[0.9,0.9,0.9],'linewidth',3,'Handlevisibility','Off')
hold on
y = squeeze(data_rr.VR_mat(ikaus,idp,2,:))';

plot(t*1000, y/max(y))
yvals = squeeze(data_rr.rSVR_mat(ikaus,idp,2,:));
tvals = squeeze(data_rr.trS_mat(ikaus,idp,2,:))*1000;
yvals = yvals/max(yvals);
plot(tvals,yvals)
%%
% code copied from calcFourierGetsignal(y',t,1,1e6,1e3,inputSignal,0,0,1)
close all
i = 1;
T_end = t(end);
Tsim = t(1:end-1);
VR = y(:,1:end-1);
freqOI = 1e6;
Window = 1e3;

figure()
plot(Tsim,VR)

Fs = 1/(Tsim(2)-Tsim(1));  % extracting sampling frequency
L = length(Tsim);           % extracting signal length
N = L;
%N = 2^nextpow2(L);             % zero padding
dF = Fs/N;
Y = fft(VR(i,:),N);          % calculating FFT from measured potential
freqs = 0:dF:Fs/2-mod(N+1,2)*dF;       % determine pos frequency compontents res F = Fs/N
f_neg=-dF:-dF:-Fs/2;   % calculating negative components
freqs = [fliplr(f_neg),freqs];  % make one array
Y = fftshift(Y);                % shift fourier transform until 0 in centre

indices = freqs>=freqOI-Window  & freqs<=freqOI+Window;
indicesneg = freqs>=-freqOI-Window  & freqs<=-freqOI+Window;
YOIpos = Y(indices);  % substract region of interest

lYOI = length(YOIpos);
figure
subplot(2,2,1);plot(freqs(indices)-freqOI,lYOI*real(YOIpos/N));
subplot(2,2,2);plot(freqs(indices)-freqOI,lYOI*imag(YOIpos/N));
subplot(2,2,3);plot(freqs(indices)-freqOI,abs(lYOI*YOIpos/N));
subplot(2,2,4);plot(freqs(indices)-freqOI,atan2(imag(YOIpos/N),real(YOIpos/N)));
title('subtraction region of interest')


%YOIneg = Y(indicesneg);
ReYOIpos = real(YOIpos);
ImYOIpos = imag(YOIpos);
% filter out unwanted signals
filterlength = min(max(ceil(length(ReYOIpos)/8),20),floor(length(ReYOIpos)/2)); % determine mean at 1Mhz from 20 or more points if size larger than number of points*2
filterArray = [1:filterlength,length(ReYOIpos)-filterlength+1:length(ReYOIpos)]; %mean static component determined from outerboundary of window of interest (contain less of components of our signal of interest)
if ~isempty(filterArray)
ReYOIpos = ReYOIpos-mean(ReYOIpos(filterArray));
ImYOIpos = ImYOIpos-mean(ImYOIpos(filterArray));
end

ReYOI = ReYOIpos;
ImYOI = ImYOIpos;
lYOI = length(ReYOI);

subplot(2,2,1);hold on;plot(freqs(indices)-freqOI,lYOI*ReYOI/N);
subplot(2,2,2);hold on;plot(freqs(indices)-freqOI,lYOI*ImYOI/N);
subplot(2,2,3);hold on;plot(freqs(indices)-freqOI,abs(lYOI*complex(ReYOI,ImYOI)/N));
subplot(2,2,4);hold on;plot(freqs(indices)-freqOI,atan2(ImYOI,ReYOI));
title('after static signal filtering')

fYOI = 2*complex(ReYOI,ImYOI)/complex(0,-1); % compensate for sin
% compensate for resampling
fYOI = length(fYOI)/N.*fYOI;

figure
subplot(2,2,1);plot(freqs(indices),real(fYOI));
subplot(2,2,2);plot(freqs(indices),imag(fYOI));
subplot(2,2,3);plot(freqs(indices),abs(fYOI));
subplot(2,2,4);plot(freqs(indices),atan2(imag(fYOI),real(fYOI)));
title('after correction for sinwave')

[S2,fout,power2] = extractSymsignal_hermit(fYOI,freqs(indices),freqOI,0,0); % make signal hermitic ==> inverse is reeel ==> no error with ifft
Sout = fftshift(fft(S2));
plot(S2)

%%
close all
clr = lines(2);
clr = clr(2,:);
clr = [0.3,0.3,0.3]+0.1;
figure()
plot(Tsim*1000,VR,'color',clr)
xlabel('time [ms]');
ylabel('\psi')
ax = gca;
ax.YAxis.Visible = 'off';
set(get(gca,'YLabel'),'visible','on')
set(gca,{'ytick','ytickLabel','box'},{[],[],'off'});
set(gcf,{'units','color','position','paperunits','papersize','Renderer'},...
    {'centimeters',[1,1,1],[1,3,3.5,3.5],'centimeters',[4,4],'Painters'})
set(findall(gcf,'type','axes'),'fontsize',9)


figure()
plot(freqs,abs(Y/N),'color',clr)
set(gca,{'xscale','yscale'},{'log','log'})
xlabel('frequency [Hz]');
ylabel('|\Psi|')
ax = gca;
ax.YAxis.Visible = 'off';
set(get(gca,'YLabel'),'visible','on')
set(gca,{'ytick','ytickLabel','box'},{[],[],'off'});
set(gca,{'xtick'},{[1,1e3,1e6]});
set(gcf,{'units','color','position','paperunits','papersize','Renderer'},...
    {'centimeters',[1,1,1],[1+5,3,3.5,3.5],'centimeters',[4,4],'Painters'})
set(findall(gcf,'type','axes'),'fontsize',9)


lYOI = length(YOIpos);
figure
subplot(2,1,1);plot(freqs(indices)-freqOI,lYOI*real(YOIpos/N),'color',clr);
subplot(2,1,2);plot(freqs(indices)-freqOI,lYOI*imag(YOIpos/N),'color',clr);
subplot(2,1,1);hold on;plot(freqs(indices)-freqOI,lYOI*ReYOI/N,':','color',min(clr+0.5,1));
ylabel('Re(\Psi)')
ax = gca;
ax.YAxis.Visible = 'off';
set(get(gca,'YLabel'),'visible','on')
set(gca,{'ytick','ytickLabel','box'},{[],[],'off'});

subplot(2,1,2);hold on;plot(freqs(indices)-freqOI,lYOI*ImYOI/N,':','color',min(clr+0.5,1));
xlabel('frequency [Hz]');
ylabel('Im(\Psi)')
ax = gca;
ax.YAxis.Visible = 'off';
set(get(gca,'YLabel'),'visible','on')
set(gca,{'ytick','ytickLabel','box'},{[],[],'off'});
set(gcf,{'units','color','position','paperunits','papersize','Renderer'},...
    {'centimeters',[1,1,1],[1+10,3,3.5,3.5],'centimeters',[4,4],'Painters'})
set(findall(gcf,'type','axes'),'fontsize',9)

figure
subplot(2,1,1);plot(freqs(indices)-freqOI,real(fYOI),'color',clr);
ylabel('Re(\Psi)')
ax = gca;
ax.YAxis.Visible = 'off';
set(get(gca,'YLabel'),'visible','on')
set(gca,{'ytick','ytickLabel','box'},{[],[],'off'});
subplot(2,1,2);plot(freqs(indices)-freqOI,imag(fYOI),'color',clr);
xlabel('frequency [Hz]');
ylabel('Im(\Psi)')
ax = gca;
ax.YAxis.Visible = 'off';
set(get(gca,'YLabel'),'visible','on')
set(gca,{'ytick','ytickLabel','box'},{[],[],'off'});
set(gcf,{'units','color','position','paperunits','papersize','Renderer'},...
    {'centimeters',[1,1,1],[1+15,3,3.5,3.5],'centimeters',[4,4],'Painters'})
set(findall(gcf,'type','axes'),'fontsize',9)

figure()
yvals = squeeze(data_rr.rSVR_mat(ikaus,idp,2,:));
tvals = squeeze(data_rr.trS_mat(ikaus,idp,2,:))*1000;
yvals = yvals/max(yvals);
Tsim2 = 0:0.1:data_rr.trS_mat(ikaus,idp,end,end)*1000;
plot(Tsim2,Ifun(Tsim2/1000),'color',[0.9,0.9,0.9],'linewidth',3,'Handlevisibility','Off')
hold on
plot(tvals,yvals,'color',clr)
xlabel('time [ms]')
ylabel("f'(t)")
ax = gca;
ax.YAxis.Visible = 'off';
set(get(gca,'YLabel'),'visible','on')
set(gca,{'ytick','ytickLabel','box'},{[],[],'off'});
set(gcf,{'units','color','position','paperunits','papersize','Renderer'},...
    {'centimeters',[1,1,1],[1+20,3,3.5,3.5],'centimeters',[4,4],'Painters'})
set(findall(gcf,'type','axes'),'fontsize',9)
%%

function [S2,fout,power2] = extractSymsignal_hermit(YOIpos,FreqsOI,freqOI,theta_us,debug_flag)
    if ~any(FreqsOI==freqOI)
        disp('caution freq of interest not in extracted frequencies ==> reconstruction?')
    else
        if mod(length(FreqsOI),2)==0
            idxfreqOI = find(FreqsOI==freqOI);
            if idxfreqOI == (length(FreqsOI)/2)+1
                posY = YOIpos(idxfreqOI+1:end);
                negY = YOIpos(2:idxfreqOI-1);
                nyquistY = YOIpos(1);
                posFreqs = FreqsOI(idxfreqOI+1:end)-freqOI;
                negFreqs = FreqsOI(2:idxfreqOI-1)-freqOI;
                nyquistFreq = FreqsOI(1)-freqOI;
                negYinterp = interp1(abs(negFreqs),negY,posFreqs,'pchip');% interpolation also filps signal => f(w) = f(-w)
                %posYinterp = interp1(posFreqs,posY,fliplr(abs(negFreqs(2:end))),'pchip'); %start from 2 contains only once nyquist
                % making Y hermitian
                % https://en.wikipedia.org/wiki/Hermitian_function ==>
                % inverse fourier is reeel
                Gmwster = conj(negYinterp);
                % because it can contain phase of sinus wave (change in
                % position of dipole wrt poi we need to correct for this
                % phase
                %-----------------------------
                %             eipsi = mean(posY./Gmwster);%e^i2theta
                %             psi = atan2(imag(eipsi),real(eipsi));%2theta
                %             corrf = complex(0,psi/2);%itheta
                %             Fmwster = exp(corrf)*Gmwster;%F*(-w)
                %             Fw = exp(-corrf)*posY;%F(w)
                %             posYnew = 1/2*(Fw+Fmwster);
                %-----------------------------
                
                modG = vecnorm(posY,2,1);
                modGster = vecnorm(Gmwster,2,1);
                mod_avg = (modG+modGster)/2;
                tppsi = atan2(imag(posY),real(posY));
                mtppsi = atan2(imag(Gmwster),real(Gmwster));
                theta_avg = theta_us;%get_thetaavg(tppsi,mtppsi)%+0*(1-mirror)*pi;
                if any(modG==0) || any(modGster==0)
                    Kw = posY.*exp(-complex(0,theta_avg));
                    Kmwster = Gmwster.*exp(complex(0,theta_avg));
                    
                else
                    Kw = mod_avg./modG.*posY.*exp(-complex(0,theta_avg));
                    Kmwster = mod_avg./modGster.*Gmwster.*exp(complex(0,theta_avg));
                    
                end
                posYnew = 1/2*(Kw+Kmwster);% 1/2*(F(w)+F*(-w))
                cYnew = YOIpos(idxfreqOI).*exp(-complex(0,theta_avg));
                nyqYnew = nyquistY.*exp(-complex(0,theta_avg));
                %-------------------------------
                %posYnew = 1/2*(posY+Gmwster);
                % cYnew = YOIpos(idxfreqOI)
                %------------------------------
                if debug_flag && PLOT
                    figure
                    subplot(2,2,1);plot(posFreqs,real(posY));
                    subplot(2,2,2);plot(posFreqs,imag(posY));
                    subplot(2,2,3);plot(posFreqs,abs(posY));
                    subplot(2,2,4);plot(posFreqs,atan2(imag(posY),real(posY)));
                    title('G(w)')
                end
                if debug_flag && PLOT
                    figure
                    subplot(2,2,1);plot(posFreqs,real(Gmwster));
                    subplot(2,2,2);plot(posFreqs,imag(Gmwster));
                    subplot(2,2,3);plot(posFreqs,abs(Gmwster));
                    subplot(2,2,4);plot(posFreqs,atan2(imag(Gmwster),real(Gmwster)));
                    title('G(w)')
                end
                negYnew = fliplr(conj(posYnew));
                Sout = [nyqYnew,negYnew,cYnew,posYnew];
                fout = [nyquistFreq ,negFreqs,0,posFreqs];
            else
                error('extract symsignal needs to be constructed still')
            end
        else
            idxfreqOI = find(FreqsOI==freqOI);
            posY = YOIpos(idxfreqOI+1:end);
            negY = YOIpos(1:idxfreqOI-1);
            posFreqs = FreqsOI(idxfreqOI+1:end)-freqOI; %G(w) = F(w)*e^-itheta
            negFreqs = FreqsOI(1:idxfreqOI-1)-freqOI;
            fout = [negFreqs,0,posFreqs];
            if any(posFreqs+fliplr(negFreqs))
                error('not symmetric???')
            end
            % making Y hermitian
            % https://en.wikipedia.org/wiki/Hermitian_function ==>
            % inverse fourier is reeel
            Gmwster = conj(fliplr(negY)); %G*(-w) = F*(-w)*e^-itheta
            % because it can contain phase of sinus wave (change in
            % position of dipole wrt poi we need to correct for this
            % phase
            %-----------------------------
%             eipsi = mean(posY./Gmwster);%e^i2theta
%             psi = atan2(imag(eipsi),real(eipsi));%2theta
%             corrf = complex(0,psi/2);%itheta
%             Fmwster = exp(corrf)*Gmwster;%F*(-w)
%             Fw = exp(-corrf)*posY;%F(w)
%             posYnew = 1/2*(Fw+Fmwster);
            %-----------------------------
            if debug_flag && PLOT
                figure
                subplot(2,2,1);plot(posFreqs,real(posY));
                subplot(2,2,2);plot(posFreqs,imag(posY));
                subplot(2,2,3);plot(posFreqs,abs(posY));
                subplot(2,2,4);plot(posFreqs,atan2(imag(posY),real(posY)));
                title('G(w)')
            end
            if debug_flag && PLOT
                figure
                subplot(2,2,1);plot(posFreqs,real(Gmwster));
                subplot(2,2,2);plot(posFreqs,imag(Gmwster));
                subplot(2,2,3);plot(posFreqs,abs(Gmwster));
                subplot(2,2,4);plot(posFreqs,atan2(imag(Gmwster),real(Gmwster)));
                title('G(w)')
            end
            modG = abs(posY);
            modGster = abs(Gmwster);
            mod_avg = (modG+modGster)/2;
            tppsi = atan2(imag(posY),real(posY));
            mtppsi = atan2(imag(Gmwster),real(Gmwster));
            theta_avg = theta_us;%get_thetaavg(tppsi,mtppsi)%+0*(1-mirror)*pi;
            if debug_flag
                diff_thetas = theta_avg-theta_us
            end
            Kw = mod_avg./modG.*posY.*exp(-complex(0,theta_avg));
            Kmwster = mod_avg./modGster.*Gmwster.*exp(complex(0,theta_avg));
            % solve nans
            if any(isnan(Kw))
                Kw(mod_avg==0&modG==0)=0;
                if any(isnan(Kw))
                    error('what to do?')
                end
            end
            if any(isnan(Kmwster))
                Kmwster(mod_avg==0&modGster==0)=0;
                if any(isnan(Kmwster))
                    error('what to do?')
                end
            end
            posYnew = 1/2*(Kw+Kmwster);% 1/2*(F(w)+F*(-w))
            cYnew = YOIpos(idxfreqOI).*exp(-complex(0,theta_avg));
            %-------------------------------
%             posYnew = 1/2*(posY+Gmwster);
%             cYnew = YOIpos(idxfreqOI)
            %------------------------------
            
            negYnew = fliplr(conj(posYnew));
            Sout = [negYnew,cYnew,posYnew];
            

        end
        if debug_flag && PLOT
            figure
            subplot(2,2,1);plot(fout,real(Sout));
            subplot(2,2,2);plot(fout,imag(Sout));
            subplot(2,2,3);plot(fout,abs(Sout));
            subplot(2,2,4);plot(fout,atan2(imag(Sout),real(Sout)));
            title('after correction for theta')
        end
    end
    ifYOI = ifft(ifftshift(Sout)); %ifft(fYOI,'symmetric');
    S2 = real(ifYOI);
    power2 = sum(abs(S2).^2);
end
