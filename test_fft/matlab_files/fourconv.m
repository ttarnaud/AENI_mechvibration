fs = 256;
flowpass=10;
N = 8;        % FFT much faster at this length
t = -N/2:1/fs:N/2-1/fs;      
h = sinc(2*t*flowpass);     % filter impulse reponse

% filter frequency response
H = fft(h);  

x = sin(2*pi*t*50) + sin(2*pi*t*1);   % input = dc (any signal will do)


Nrep = 100;      % number of trials to average

t0 = clock;      % latch the current time
for i=1:Nrep, y1 = conv(x,h); end      % Direct convolution
%for i=1:Nrep, y2 = ifft(fft(x) .* H); end % FFT convolution
t1 = etime(clock,t0)*1000; % elapsed time in msec


disp(sprintf('Average time = %0.2f msec\n',t1/Nrep));

%%
figure,
subplot(3,1,1);plot(x);
subplot(3,1,2);plot(y1);
subplot(3,1,3);plot(y2);