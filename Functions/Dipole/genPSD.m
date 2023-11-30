function [PSD,f] = genPSD(signal,Tsim,Fs,type)
N = length(Tsim);
%important that signal is limited to fixed period: if Tsim 0:0.1:1 Fs =
%10Hz and number of samples is 11 if signal is sin(2pi t) for this set the
%value of 0 and 1 will be the same. If we use the fft algortihm an artefact
%will be created because its is based on an infinte signal assumption.
%however in current form the signal is not perfect sinusoidal if extended
%infinitly . start each period will contain 0 twice
%therefore crop signal length = artefact reduction
if abs(signal(end)-signal(1))<max(signal)*1e-9   
signal = signal(1:end-1);
Tsim = Tsim(1:end-1);
end
N = length(Tsim);           % extracting signal length
dF = Fs/N;
%N = 2^nextpow2(L);             % zero padding  can be used to increase
%resolution of frequency domain however at the cost of a sinc artefact
%(zeropadding = muliplication of squarewave to time signal or is equal to
%convulotion of sinc with fourier transform of input signal)
Y = fft(signal,N);          % calculating FFT from measured potential 
f = 0:dF:Fs/2;             % determine pos frequency compontents res F = Fs/N when N is even the Nyquist frequencie is included once when N uneven Nyquist not included
yright = Y(1:floor(N/2)+1); %if even number of sample points only one FS/2 value inlcuded this is at N/2+1 (part of negative freq domain), when uneven number Fs/2 (near Fs/2 depending of dF) is both in pos and neg domain
PSD = 1/N*abs(yright).^2; %Normalization not squared see: https://math.stackexchange.com/questions/346894/prove-of-the-parsevals-theorem-for-discrete-fourier-transform-dft
if mod(N,2)==0
    PSD(2:end-1) = 2*PSD(2:end-1);
else
    PSD(2:end) = 2*PSD(2:end);
end

% good source with explanation on normalization https://dsp.stackexchange.com/questions/11376/why-are-magnitudes-normalised-during-synthesis-idft-not-analysis-dft

switch lower(type)
    case'db'
        PSD = 10*log10(PSD);
    case'db/freq'
        PSD = 10*log10(PSD/Fs);
    otherwise
        PSD = PSD;
end

end