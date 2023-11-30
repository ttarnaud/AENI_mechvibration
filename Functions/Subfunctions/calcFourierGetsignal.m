function [Out,Power,time] =  calcFourierGetsignal(VR,Tsim,POI,freqOI,Window,inputSignal,mirror,theta_US,PLOT,varargin)
infostr = {'Single sided',['Single sided', 'elongated signal'],'2 sided', 'Real 2 sided','Imag 2 sided','Phase 2 sided'};
Titles = {'Complete spectrum','0-1kHz',['around',num2str(freqOI*10^-6),' MHz']};
Fignr = get(gcf,'number');
method = 'compare';
inputPLOT = PLOT;
debug_flag =0;
PtP = 1;
% change accepted: this change makes that the last point is not the same as
% the first (this gave errors in fft as it disrupted the periodicity)
% see genPSD. Because signal in investbiologicalNoise contains both 0 and
% Tend and we are working a pure sine wave is used: crop is necessary
T_end = Tsim(end);
Tsim = Tsim(1:end-1);
VR = VR(:,1:end-1);

    
for i=POI

if iscell(varargin) && length(varargin) == 1 && iscell(varargin{1})
    varargin = varargin{1};
end
if mod(length(varargin)+1,2)
    if ~isempty(varargin)
        if any(strcmpi(varargin,'POIstoPLOT'))
            PtP = varargin{find(strcmpi(varargin,'POIstoPLOT'))+1};
        end
    end
end
if any(PtP==i) && inputPLOT
    PLOT=1;
else
    PLOT=0;
end
    
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
if debug_flag && PLOT
    figure
    subplot(2,4,1);plot(freqs,real(Y/N));xlim([-1000,1000])
    subplot(2,4,2);plot(freqs,imag(Y/N));xlim([-1000,1000])
    subplot(2,4,5);plot(freqs,abs(Y/N));xlim([-1000,1000])
    subplot(2,4,6);plot(freqs,atan2(imag(Y/N),real(Y/N)));xlim([-1000,1000])
    subplot(2,4,3);plot(freqs,real(Y/N));xlim([-1000,1000]+freqOI)
    subplot(2,4,4);plot(freqs,imag(Y/N));xlim([-1000,1000]+freqOI)
    subplot(2,4,7);plot(freqs,abs(Y/N));xlim([-1000,1000]+freqOI)
    subplot(2,4,8);plot(freqs,atan2(imag(Y/N),real(Y/N)));xlim([-1000,1000]+freqOI)
    title('before any adjustments (only Tsim(1:end-1))')
    mtit('before any adjustments (only Tsim(1:end-1))')
end
    
indices = freqs>=freqOI-Window  & freqs<=freqOI+Window;
indicesneg = freqs>=-freqOI-Window  & freqs<=-freqOI+Window;
YOIpos = Y(indices);  % substract region of interest

%signal power
%calc power https://en.wikipedia.org/wiki/Parseval%27s_theorem power should
%be equal to power2;
Power(i) = 1/length(YOIpos).*sum(abs(YOIpos).^2);
if debug_flag && PLOT
    figure
    subplot(2,2,1);plot(freqs(indices),real(YOIpos/N));
    subplot(2,2,2);plot(freqs(indices),imag(YOIpos/N));
    subplot(2,2,3);plot(freqs(indices),abs(YOIpos/N));
    subplot(2,2,4);plot(freqs(indices),atan2(imag(YOIpos/N),real(YOIpos/N)));
    title('subtraction region of interest')
end
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

if debug_flag && PLOT
    figure
    subplot(2,2,1);plot(freqs(indices),ReYOI);
    subplot(2,2,2);plot(freqs(indices),ImYOI);
    subplot(2,2,3);plot(freqs(indices),abs(complex(ReYOI,ImYOI)));
    subplot(2,2,4);plot(freqs(indices),atan2(ImYOI,ReYOI));
    title('after static signal filtering')
end
if mod(length(ReYOIpos),2)
ReYOI = ReYOI(1:end-1);  %make even => contains one time nyquist
ImYOI = ImYOI(1:end-1);
indices = find(indices);
indices = indices(1:end-1);
end

%dt = 1/(2*Window);
dt = 1/(dF*length(ReYOI));
%time{i} = 0:dt:T_end;
time{i} = (0:length(ReYOI)-1).*dt; 
if ~(time{i}(end)==T_end||time{i}(end)==T_end-dt)
    error('check time axis reconstructed signal')
end
fYOI = 2*complex(ReYOI,ImYOI)/complex(0,-1); % compensate for sin
% compensate for resampling
fYOI = length(fYOI)/N.*fYOI;
%calc power https://en.wikipedia.org/wiki/Parseval%27s_theorem power should
%be equal to power2;
power(i) = 1/length(fYOI).*sum(abs(fYOI).^2);
if debug_flag && PLOT
    figure
    subplot(2,2,1);plot(freqs(indices),real(fYOI));
    subplot(2,2,2);plot(freqs(indices),imag(fYOI));
    subplot(2,2,3);plot(freqs(indices),abs(fYOI));
    subplot(2,2,4);plot(freqs(indices),atan2(imag(fYOI),real(fYOI)));
    title('after correction for sinwave')
end
if any(fYOI)
    switch lower(method)
        case 'hermitian'
            [S2,fout,power2] = extractSymsignal_hermit(fYOI,freqs(indices),freqOI,theta_US,debug_flag); % make signal hermitic ==> inverse is reeel ==> no error with ifft
        case 'real'
            [S2,fout,power2] = extractSymsignal_real(fYOI,freqs(indices),freqOI,power(i),theta_US,debug_flag);
            
        case 'compare'
            [S2,fout,power2] = extractSymsignal_hermit(fYOI,freqs(indices),freqOI,theta_US,debug_flag); % make signal hermitic ==> inverse is reeel ==> no error with ifft
            [S2_re_sc,~,power2,S2_re] = extractSymsignal_real(fYOI,freqs(indices),freqOI,power(i),theta_US,debug_flag);
            
    end
else
    %if frequency compenents all zero error in extract functions above.
    S2 = fYOI; fout = freqs(indices);
    S2_re_sc = fYOI; S2_re = fYOI;
end
    
Out{i} = S2;
Sout = fftshift(fft(S2));



if PLOT
    if debug_flag
        Fignr = get(gcf,'number');
    end
figure(Fignr+i)
subplot(2,2,1)
plot(fout,real(Sout),'b')
if strcmpi(method,'compare')
    hold on
    plot(fout,real(fftshift(fft(S2_re))),'r--','DisplayName','RealRecon')
    plot(fout,real(fftshift(fft(S2_re_sc))),'g--','DisplayName','RealRecon+scale')
    hold off
end
ylabel('Re(Fx)')
xlabel('f [Hz]')
subplot(2,2,2)
plot(fout,imag(Sout),'b')
if strcmpi(method,'compare')
    hold on
    plot(fout,imag(fftshift(fft(S2_re))),'r--','DisplayName','RealRecon')
    plot(fout,imag(fftshift(fft(S2_re_sc))),'g--','DisplayName','RealRecon+scale')
    hold off
end
ylabel('Im(Fx)')
xlabel('f [Hz]')
if isempty(inputSignal)
subplot(2,2,[3,4])
plot(time{i}*1e3,S2,'b')
if strcmpi(method,'compare')
    hold on
    plot(time{i}*1e3,S2_re,'r--','DisplayName','RealRecon')
    plot(time{i}*1e3,S2_re_sc,'g--','DisplayName','RealRecon+scale')
    hold off
end
ylabel('V')
xlabel('time [ms]')
else
    subplot(2,2,3)
    plot(time{i}*1e3,S2,'b')
    if strcmpi(method,'compare')
    hold on
    plot(time{i}*1e3,S2_re,'r--','DisplayName','RealRecon')
    plot(time{i}*1e3,S2_re_sc,'g--','DisplayName','RealRecon+scale')
    hold off
    end
    ylabel('V')
    xlabel('Time [ms]') 
    subplot(2,2,4)
    plot(time{i}*1e3,(-1)^mirror(i)*S2./max(abs(S2)),'b','DisplayName','Reconstructed signal')
    hold on
    plot(inputSignal(1,:)*1e3,inputSignal(2,:)./max(abs(inputSignal(2,:))),'k','DisplayName','Input signal')
    if strcmpi(method,'compare')
    
    plot(time{i}*1e3,(-1)^mirror(i).*S2_re./max(abs(S2_re)),'r--','linewidth',2,'DisplayName','RealRecon')
    plot(time{i}*1e3,(-1)^mirror(i).*S2_re_sc./max(abs(S2_re_sc)),'g--','DisplayName','RealRecon+scale')
    
    end
    hold off
    legend('show')
    ylabel('normalized signal')
    xlabel('Time [ms]')
end
mtit(['POI ',num2str(i)])
if i==1
Fignr2 = Fignr+max(PtP)+1;
end
if PLOT
figure(Fignr2)
subplot(length(PtP),2,1+2*(i-1))
plot(time{i}*1e3,S2*1e6)
ylabel({['\bf POI ',num2str(i)],'','V [pV]'})
xlabel('Time [ms]')
subplot(length(PtP),2,2+2*(i-1))
plot(time{i}*1e3,(-1)^mirror(i).*S2./max(abs(S2)),'DisplayName','Reconstructed signal')
hold on
plot(inputSignal(1,:)*1e3,inputSignal(2,:)./max(abs(inputSignal(2,:))),'DisplayName','Input signal')
hold off
legend('show')
ylabel('normalized signal')
xlabel('Time [ms]')
pause(0.01)
end
end
end
Fignr = get(gcf,'number');
if PLOT && debug_flag
for ifig =PtP
    XLIMS = [0,1e3;(freqOI-Window)*1e-6,(freqOI+Window)*1e-6];
    figure(Fignr+ifig)
    for j = 1:6
        for isub = 1:3
            subplot(6,3,(j-1)*3+isub);
            nr = ifig;
            if j == 2
                N = 2^nextpow2(L); % to get higher resolution willl giv sinc function
            else
                N = L;
            end
            Y = fft(VR(nr,:),N);
            
            pfreqs = 0:dF:Fs/2;
            if j>=3
                f_neg=-dF:-dF:-Fs/2;
                f_pos = 0:dF:Fs/2-mod(N+1,2)*dF;
                pfreqs = [fliplr(f_neg),f_pos];
                Y = fftshift(Y);
                XLIMS = [-Window,Window;XLIMS(2,:)];
            end
            if isub == 3
                pfreqs = pfreqs/1e6;
            end
            if j==4
                plot(pfreqs,real(Y(1:length(pfreqs))))
            elseif j==5
                plot(pfreqs,imag(Y(1:length(pfreqs))))
            elseif j==6
                plot(pfreqs,angle(Y(1:length(pfreqs))))
            else
                plot(pfreqs,abs(Y(1:length(pfreqs))))
            end
          
            xlabel('f [Hz]')
            if j == 1
                title(Titles{isub})
            end
            if isub == 1
                ylabel({['\bf', infostr{j}],'','Fx'})
            end
            if isub>1
                xlim(XLIMS(isub-1,:))
            end
            if isub == 3
                %set(gca,'XTickLabel',arrayfun(@num2str,get(gca,'xTick')*1e-6,'UniformOutput',false))
                xlabel('f [MHz]')
            end            
        end       
        
    end
    mtit(['POI',num2str(nr)])
    pause(0.01)
end
end
function [S2_sc,fout,power,S2] = extractSymsignal_real(fYOI,FreqsOI,freqOI,power_init,theta_us,debug_flag)
        idxfreqOI = find(FreqsOI==freqOI);
        posFreqs = FreqsOI(idxfreqOI+1:end)-freqOI;
        negFreqs = FreqsOI(1:idxfreqOI-1)-freqOI;
        fout = [negFreqs,0,posFreqs];
        fYOI = fYOI.*exp(-complex(0,theta_us));
        if debug_flag && PLOT
            figure
            subplot(2,2,1);plot(fout,real(fYOI));
            subplot(2,2,2);plot(fout,imag(fYOI));
            subplot(2,2,3);plot(fout,abs(fYOI));
            subplot(2,2,4);plot(fout,atan2(imag(fYOI),real(fYOI)));
            title('after correction for theta sinwave')
        end
        
        ifYOI = ifft(ifftshift(fYOI)); %ifft(fYOI,'symmetric');
        S2 = real(ifYOI);
        power = sum(abs(S2).^2);
        ratioPower = power_init/power;
        S2_sc = sqrt(ratioPower).*S2;
        if debug_flag && PLOT
            fYOI = fftshift(fft(S2_sc));
            figure
            subplot(2,2,1);plot(fout,real(fYOI));
            subplot(2,2,2);plot(fout,imag(fYOI));
            subplot(2,2,3);plot(fout,abs(fYOI));
            subplot(2,2,4);plot(fout,atan2(imag(fYOI),real(fYOI)));
            title('after correction power in real reconstreuction')
        end
end
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
    function theta_avg = get_thetaavg(tppsi,mtppsi)
        %https://en.wikipedia.org/wiki/Mean_of_circular_quantities 
        normunitvec = vecnorm((exp(complex(0,tppsi))+exp(complex(0,mtppsi)))/2,2,1);
        warn_flag = any(normunitvec<0.5);
        if warn_flag & 0
            disp(['norm average unit vector smaller than 0.5 on position',num2str(find(normunitvec<0.5))]);
        end
        mre = (cos(tppsi)+cos(mtppsi))/2;
        mim = (sin(tppsi)+sin(mtppsi))/2;
        psipnpi = atan2(mim,mre);
        theta = tppsi-psipnpi;
        theta_re = cos(theta);
        theta_im = sin(theta);
        normunitvec = vecnorm(mean(exp(complex(0,theta))),2,1);
        warn_flag2 = any(normunitvec<0.5);
        if warn_flag2 &0
            disp(['norm average unit vector smaller than 0.5 on position',num2str(find(normunitvec<0.5))]);
        end
        theta_avg = atan2(mean(theta_im),mean(theta_re));
    end
end