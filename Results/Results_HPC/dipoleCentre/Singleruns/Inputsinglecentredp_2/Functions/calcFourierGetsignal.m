function [Out,power,time] =  calcFourierGetsignal(VR,Tsim,POI,freqOI,Window,inputSignal,mirror,PLOT)
infostr = {'Single sided',['Single sided', 'elongated signal'],'2 sided', 'Real 2 sided','Imag 2 sided','Phase 2 sided'};
Titles = {'Complete spectrum','0-1kHz',['around',num2str(freqOI*10^-6),' MHz']};
Fignr = get(gcf,'number');
for i=POI
Fs = 1/(Tsim(2)-Tsim(1));  % extracting sampling frequency
L = length(Tsim);           % extracting signal length
N = L;
%N = 2^nextpow2(L);             % zero padding
Y = fft(VR(i,:),N);          % calculating FFT from measured potential
dF = 1/Tsim(end);
freqs = 0:dF:Fs/2;       % determine pos frequency compontents res F = Fs/N
f_neg=-dF:-dF:-Fs/2+mod(N+1,2)*dF;   % calculating negative components
freqs = [fliplr(f_neg),freqs];  % make one array
Y = fftshift(Y);                % shift fourier transform until 0 in centre
indices = freqs>=freqOI-Window  & freqs<=freqOI+Window;
indicesneg = freqs>=-freqOI-Window  & freqs<=-freqOI+Window;
YOIpos = Y(indices);  % substract region of interest
%YOIneg = Y(indicesneg);
ReYOIpos = real(YOIpos);
ImYOIpos = imag(YOIpos);
% filter out unwanted signals
ReYOIpos = ReYOIpos-mean(ReYOIpos);
ImYOIpos = ImYOIpos-mean(ImYOIpos);
ReYOI = ReYOIpos;
ImYOI = ImYOIpos;



dt = 1/(2*Window);
time{i} = 0:dt:Tsim(end);
fYOI = complex(ReYOI,ImYOI)/complex(0,-1);
[Sout,fout] = extractSymsignal(fYOI,freqs(indices),freqOI); % make signal symmetric ==> no error with ifft
fYOI = ifftshift(Sout);
ifYOI = ifft(fYOI); %ifft(fYOI,'symmetric');
S2 = 2*length(YOIpos)/N*real(ifYOI);
Out{i} = S2;

%calc power
power(i) = sum(1/N*abs(fYOI).^2);


if PLOT
figure(Fignr+i)
subplot(2,2,1)
plot(fout,real(Sout))
ylabel('Re(Fx)')
xlabel('f [Hz]')
subplot(2,2,2)
plot(fout,imag(Sout))
ylabel('Im(Fx)')
xlabel('f [Hz]')
if isempty(inputSignal)
subplot(2,2,[3,4])
plot(time{i}*1e3,S2)
ylabel('V')
xlabel('time [ms]')
else
    subplot(2,2,3)
    plot(time{i}*1e3,S2)
    ylabel('V')
    xlabel('Time [ms]') 
    subplot(2,2,4)
    plot(time{i}*1e3,(-1)^mirror(i)*S2./max(abs(S2)),'DisplayName','Reconstructed signal')
    hold on
    plot(inputSignal(1,:)*1e3,inputSignal(2,:)./max(abs(inputSignal(2,:))),'DisplayName','Input signal')
    hold off
    legend('show')
    ylabel('normalized signal')
    xlabel('Time [ms]')
end
mtit(['POI ',num2str(i)])

figure(Fignr+max(POI)+1)
subplot(3,2,1+2*(i-1))
plot(time{i}*1e3,S2*1e6)
ylabel({['\bf POI ',num2str(i)],'','V [pV]'})
xlabel('Time [ms]')
subplot(3,2,2+2*(i-1))
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

if PLOT
for ifig =POI
    XLIMS = [0,1e3;(freqOI-Window)*1e-6,(freqOI+Window)*1e-6];
    figure(Fignr+max(POI)+1+ifig)
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
                f_neg=-dF:-dF:-Fs/2+mod(N+1,2)*dF;
                pfreqs = [fliplr(f_neg),pfreqs];
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
function [Sout,fout] = extractSymsignal(YOIpos,FreqsOI,freqOI)
    if ~any(FreqsOI==freqOI)
        disp('caution freq of interest not in extracted frequencies ==> reconstruction?')
    else
        if mod(length(FreqsOI),2)==0
            idxfreqOI = find(FreqsOI==freqOI);
            if idxfreqOI == floor(length(FreqsOI)/2)
                posY = YOIpos(idxfreqOI+1:end);
                negY = YOIpos(1:idxfreqOI-1);
                posFreqs = FreqsOI(idxfreqOI+1:end)-freqOI;
                negFreqs = FreqsOI(1:idxfreqOI-1)-freqOI;
                negYinterp = interp1(abs(negFreqs),negY,posFreqs,'pchip');
                posYnew = 1/2*(posY+conj(negYinterp));
                negYnew = fliplr(conj(posYnew));
                Sout = [negYnew,YOIpos(idxfreqOI),posYnew];
                fout = [-1*fliplr(posFreqs),0,posFreqs];
            else
                error('extract symsignal needs to be constructed still')
            end
        else
            idxfreqOI = find(FreqsOI==freqOI);
            posY = YOIpos(idxfreqOI+1:end);
            negY = YOIpos(1:idxfreqOI-1);
            posFreqs = FreqsOI(idxfreqOI+1:end)-freqOI;
            negFreqs = FreqsOI(1:idxfreqOI-1)-freqOI;
            if any(posFreqs+fliplr(negFreqs))
                error('not symmetric???')
            end
            posYnew = 1/2*(posY+conj(fliplr(negY)));
            negYnew = fliplr(conj(posYnew));
            Sout = [negYnew,YOIpos(idxfreqOI),posYnew];
            fout = [negFreqs,0,posFreqs];
        end
    end
end
end