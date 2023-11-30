function [Iarray,extra_opt] = getIarray(type,Itype,Tsim,N,d,I,sigmaI,rselect)
%this function creates an array 'Iarray' with same lengths of number of dipoles N
%that contains the current fluctuation over time Tsim
%Inputs:
%   type =
%   'randomdelay','randomi','randominewgen','randomirandomdelay','sync',
%   determines how signals vary over time
%   Itype = 'same','normal','lognormal'. determines how maximum current is
%   distributed over space
%   Tsim: duration of current signal
%   N: number of dipoles this run
%   d: distance current source an sink
%   sigmaI: standard deviation of max I
%   I: max current during this time signal
%   rselect: currenttraces to select from loaded input (in randomi and randomdelay)
%   Itime: function

Display = 0; %display progess
Pflag = 0;   %progress only certain percentages
Iarray = zeros(N,length(Tsim));
PC_flag = 0 && ~Display; %parallelcomputing
Itime_type = 'hipPyr25us1';
debug_flag = 0;
extra_opt = struct();

switch lower(type)
    case lower('artificialPSD')
        % the currents are derived from a aritificial PSD  based on
        % observed trends in literature for low frequency(<100Hz) and model
        % ouputs (Neuron) High frequency. A plausible sprectrum is created.
        % This way minimizes interpolation errors. see An_spec script for
        % more details
        L = length(Tsim);
        Fs = 1/(Tsim(2)-Tsim(1));
        dF = Fs/L;
        posf = 0:dF:Fs/2; 
        
        wn = 14;
        theta = 0.7*rand()+0.3;         %max is reached between wn and theta*wn
        w2 = 10^3;                      %breakpoint higher powerlaw
        alpha1 = rand()+2;
        alpha2 = 2*rand()+2;
        %noise2 =@(x) 10.^(randn(1,length(x)));  %different noise methods
        %were tested see script An_spec
        noise3 = @(x) lognrnd(0,0.5,1,length(x));
        y = @(x) (1+x/(theta*wn))./((1+(x/wn).^alpha1).*(1+(x/w2).^alpha2));
        fun = @(x) y(x).*noise3(x); %this function determines artificial PSD right hand sided
        
        %store parameters
        extra_opt.alpha1 = alpha1; extra_opt.alpha2 = alpha2; extra_opt.theta = theta;
        if mod(L,2)==0
            
            PSD = [fun(posf(1)),fun(posf(2:end-1))/2,...
                fun(posf(end)),fun(fliplr(posf(2:end-1)))/2]; %same way as components would be stored by using fft
            if debug_flag
            abs_freqs=[posf,fliplr(posf(2:end-1))];
            PSDpos = PSD(1:length(posf));
            PSDpos(2:end-1) = 2*PSDpos(2:end-1);
            end
            
        else
            
            PSD = [fun(posf(1)),fun(posf(2:end))/2,...
                fun(fliplr(posf(2:end)))/2];
            if debug_flag
            PSDpos = PSD(1:length(posf));
            PSDpos(2:end) = 2*PSDpos(2:end);
            abs_freqs=[posf,fliplr(posf(2:end))];
            end
        end
        
        s0dft = sqrt(L*PSD); %PSD to magnitude of rotor
        
        Phifk = 2*pi*rand(N,L);   %random phases
        Zfk = s0dft.*exp(complex(0,Phifk)); %add random phases while same frequency content maintained
        Iarray = ifft(Zfk,[],2,'symmetric'); %symmetric to get real output however also some artefacts but not important here (see calcgetfourier signal)
        Iarray = Iarray./max(abs(Iarray),[],2);
        
        if debug_flag
            %create plots in case debug needed
            for iN=1:N
                
                [N_p4(iN,:),N_f4(iN,:)] = genPSD(Iarray(iN,:),Tsim,Fs,'db');
            end         
            iNs = ranperm(N,min(10,N));
            figure
            subplot(2,2,[1,2])
            for iN=iNs
                plot(Tsim,Iarray(iN,:),'displayName','Iarray_temp')
                hold on
            end
            legend('show')
            ylim([-1.5,1.5])
            subplot(2,2,3)
            plot(posf,10*log10(PSDpos),'displayName','fft:Mothersignal')
            
            set(gca,'xscale','log')
            legend('show')
            subplot(2,2,4)
            for iN=iNs
                plot(N_f4(iN,:),N_p4(iN,:),'displayName','Mothersignal+random phases')
                hold on
            end
            hold off
            
        end
    case lower('artificialPSD_fixed')
        % the currents are derived from a aritificial PSD  based on
        % observed trends in literature for low frequency(<100Hz) and model
        % ouputs (Neuron) High frequency. A plausible sprectrum is created.
        % This way minimizes interpolation errors. see An_spec script for
        % more details. the powerlaws are fixed 
        L = length(Tsim);
        Fs = 1/(Tsim(2)-Tsim(1));
        dF = Fs/L;
        posf = 0:dF:Fs/2;
        
        wn = 14;
        theta = 0.6;         %max is reached between wn and theta*wn
        w2 = 10^3;                      %breakpoint higher powerlaw
        alpha1 = 2;
        alpha2 = 2;
        %noise2 =@(x) 10.^(randn(1,length(x)));  %different noise methods
        %were tested see script An_spec
        noise3 = @(x) lognrnd(0,0.5,1,length(x));
        y = @(x) (1+x/(theta*wn))./((1+(x/wn).^alpha1).*(1+(x/w2).^alpha2));
        fun = @(x) y(x).*noise3(x); %this function determines artificial PSD right hand sided
        
        %store parameters
        extra_opt.alpha1 = alpha1; extra_opt.alpha2 = alpha2; extra_opt.theta = theta;
        if mod(L,2)==0
            
            PSD = [fun(posf(1)),fun(posf(2:end-1))/2,...
                fun(posf(end)),fun(fliplr(posf(2:end-1)))/2]; %same way as components would be stored by using fft
            if debug_flag
                abs_freqs=[posf,fliplr(posf(2:end-1))];
                PSDpos = PSD(1:length(posf));
                PSDpos(2:end-1) = 2*PSDpos(2:end-1);
            end
            
        else
            
            PSD = [fun(posf(1)),fun(posf(2:end))/2,...
                fun(fliplr(posf(2:end)))/2];
            if debug_flag
                PSDpos = PSD(1:length(posf));
                PSDpos(2:end) = 2*PSDpos(2:end);
                abs_freqs=[posf,fliplr(posf(2:end))];
            end
        end
        
        s0dft = sqrt(L*PSD); %PSD to magnitude of rotor
        
        Phifk = 2*pi*rand(N,L);   %random phases
        Zfk = s0dft.*exp(complex(0,Phifk)); %add random phases while same frequency content maintained
        Iarray = ifft(Zfk,[],2,'symmetric'); %symmetric to get real output however also some artefacts but not important here (see calcgetfourier signal)
        Iarray = Iarray./max(abs(Iarray),[],2);
        
        if debug_flag
            %create plots in case debug needed
            for iN=1:N
                
                [N_p4(iN,:),N_f4(iN,:)] = genPSD(Iarray(iN,:),Tsim,Fs,'db');
            end
            iNs = ranperm(N,min(10,N));
            figure
            subplot(2,2,[1,2])
            for iN=iNs
                plot(Tsim,Iarray(iN,:),'displayName','Iarray_temp')
                hold on
            end
            legend('show')
            ylim([-1.5,1.5])
            subplot(2,2,3)
            plot(posf,10*log10(PSDpos),'displayName','fft:Mothersignal')
            
            set(gca,'xscale','log')
            legend('show')
            subplot(2,2,4)
            for iN=iNs
                plot(N_f4(iN,:),N_p4(iN,:),'displayName','Mothersignal+random phases')
                hold on
            end
            hold off
            
        end
    case lower('artificialPSD_fixed2')
        % the currents are derived from a aritificial PSD  based on
        % observed trends in literature for low frequency(<100Hz) and model
        % ouputs (Neuron) High frequency. A plausible sprectrum is created.
        % This way minimizes interpolation errors. see An_spec script for
        % more details. the powerlaws are fixed
        L = length(Tsim);
        Fs = 1/(Tsim(2)-Tsim(1));
        dF = Fs/L;
        posf = 0:dF:Fs/2;
        
        wn = 14;
        theta = 0.6;         %max is reached between wn and theta*wn
        w2 = 10^3;                      %breakpoint higher powerlaw
        alpha1 = 2;
        alpha2 = 3;
        %noise2 =@(x) 10.^(randn(1,length(x)));  %different noise methods
        %were tested see script An_spec
        noise3 = @(x) lognrnd(0,0.5,1,length(x));
        y = @(x) (1+x/(theta*wn))./((1+(x/wn).^alpha1).*(1+(x/w2).^alpha2));
        fun = @(x) y(x).*noise3(x); %this function determines artificial PSD right hand sided
        
        %store parameters
        extra_opt.alpha1 = alpha1; extra_opt.alpha2 = alpha2; extra_opt.theta = theta;
        if mod(L,2)==0
            
            PSD = [fun(posf(1)),fun(posf(2:end-1))/2,...
                fun(posf(end)),fun(fliplr(posf(2:end-1)))/2]; %same way as components would be stored by using fft
            if debug_flag
                abs_freqs=[posf,fliplr(posf(2:end-1))];
                PSDpos = PSD(1:length(posf));
                PSDpos(2:end-1) = 2*PSDpos(2:end-1);
            end
            
        else
            
            PSD = [fun(posf(1)),fun(posf(2:end))/2,...
                fun(fliplr(posf(2:end)))/2];
            if debug_flag
                PSDpos = PSD(1:length(posf));
                PSDpos(2:end) = 2*PSDpos(2:end);
                abs_freqs=[posf,fliplr(posf(2:end))];
            end
        end
        
        s0dft = sqrt(L*PSD); %PSD to magnitude of rotor
        
        Phifk = 2*pi*rand(N,L);   %random phases
        Zfk = s0dft.*exp(complex(0,Phifk)); %add random phases while same frequency content maintained
        Iarray = ifft(Zfk,[],2,'symmetric'); %symmetric to get real output however also some artefacts but not important here (see calcgetfourier signal)
        Iarray = Iarray./max(abs(Iarray),[],2);
        
        if debug_flag
            %create plots in case debug needed
            for iN=1:N
                
                [N_p4(iN,:),N_f4(iN,:)] = genPSD(Iarray(iN,:),Tsim,Fs,'db');
            end
            iNs = ranperm(N,min(10,N));
            figure
            subplot(2,2,[1,2])
            for iN=iNs
                plot(Tsim,Iarray(iN,:),'displayName','Iarray_temp')
                hold on
            end
            legend('show')
            ylim([-1.5,1.5])
            subplot(2,2,3)
            plot(posf,10*log10(PSDpos),'displayName','fft:Mothersignal')
            
            set(gca,'xscale','log')
            legend('show')
            subplot(2,2,4)
            for iN=iNs
                plot(N_f4(iN,:),N_p4(iN,:),'displayName','Mothersignal+random phases')
                hold on
            end
            hold off
            
        end
    case lower('artificialPSD_fixed3')
        % the currents are derived from a aritificial PSD  based on
        % observed trends in literature for low frequency(<100Hz) and model
        % ouputs (Neuron) High frequency. A plausible sprectrum is created.
        % This way minimizes interpolation errors. see An_spec script for
        % more details. the powerlaws are fixed
        L = length(Tsim);
        Fs = 1/(Tsim(2)-Tsim(1));
        dF = Fs/L;
        posf = 0:dF:Fs/2;
        
        wn = 14;
        theta = 0.6;         %max is reached between wn and theta*wn
        w2 = 10^3;                      %breakpoint higher powerlaw
        alpha1 = 2;
        alpha2 = 4;
        %noise2 =@(x) 10.^(randn(1,length(x)));  %different noise methods
        %were tested see script An_spec
        noise3 = @(x) lognrnd(0,0.5,1,length(x));
        y = @(x) (1+x/(theta*wn))./((1+(x/wn).^alpha1).*(1+(x/w2).^alpha2));
        fun = @(x) y(x).*noise3(x); %this function determines artificial PSD right hand sided
        
        %store parameters
        extra_opt.alpha1 = alpha1; extra_opt.alpha2 = alpha2; extra_opt.theta = theta;
        if mod(L,2)==0
            
            PSD = [fun(posf(1)),fun(posf(2:end-1))/2,...
                fun(posf(end)),fun(fliplr(posf(2:end-1)))/2]; %same way as components would be stored by using fft
            if debug_flag
                abs_freqs=[posf,fliplr(posf(2:end-1))];
                PSDpos = PSD(1:length(posf));
                PSDpos(2:end-1) = 2*PSDpos(2:end-1);
            end
            
        else
            
            PSD = [fun(posf(1)),fun(posf(2:end))/2,...
                fun(fliplr(posf(2:end)))/2];
            if debug_flag
                PSDpos = PSD(1:length(posf));
                PSDpos(2:end) = 2*PSDpos(2:end);
                abs_freqs=[posf,fliplr(posf(2:end))];
            end
        end
        
        s0dft = sqrt(L*PSD); %PSD to magnitude of rotor
        
        Phifk = 2*pi*rand(N,L);   %random phases
        Zfk = s0dft.*exp(complex(0,Phifk)); %add random phases while same frequency content maintained
        Iarray = ifft(Zfk,[],2,'symmetric'); %symmetric to get real output however also some artefacts but not important here (see calcgetfourier signal)
        Iarray = Iarray./max(abs(Iarray),[],2);
        
        if debug_flag
            %create plots in case debug needed
            for iN=1:N
                
                [N_p4(iN,:),N_f4(iN,:)] = genPSD(Iarray(iN,:),Tsim,Fs,'db');
            end
            iNs = ranperm(N,min(10,N));
            figure
            subplot(2,2,[1,2])
            for iN=iNs
                plot(Tsim,Iarray(iN,:),'displayName','Iarray_temp')
                hold on
            end
            legend('show')
            ylim([-1.5,1.5])
            subplot(2,2,3)
            plot(posf,10*log10(PSDpos),'displayName','fft:Mothersignal')
            
            set(gca,'xscale','log')
            legend('show')
            subplot(2,2,4)
            for iN=iNs
                plot(N_f4(iN,:),N_p4(iN,:),'displayName','Mothersignal+random phases')
                hold on
            end
            hold off
            
        end
    case lower('artificialPSD_fixed4')
        % the currents are derived from a aritificial PSD  based on
        % observed trends in literature for low frequency(<100Hz) and model
        % ouputs (Neuron) High frequency. A plausible sprectrum is created.
        % This way minimizes interpolation errors. see An_spec script for
        % more details. the powerlaws are fixed
        L = length(Tsim);
        Fs = 1/(Tsim(2)-Tsim(1));
        dF = Fs/L;
        posf = 0:dF:Fs/2;
        
        wn = 14;
        theta = 0.6;         %max is reached between wn and theta*wn
        w2 = 10^3;                      %breakpoint higher powerlaw
        alpha1 = 2;
        alpha2 = 5;
        %noise2 =@(x) 10.^(randn(1,length(x)));  %different noise methods
        %were tested see script An_spec
        noise3 = @(x) lognrnd(0,0.5,1,length(x));
        y = @(x) (1+x/(theta*wn))./((1+(x/wn).^alpha1).*(1+(x/w2).^alpha2));
        fun = @(x) y(x).*noise3(x); %this function determines artificial PSD right hand sided
        
        %store parameters
        extra_opt.alpha1 = alpha1; extra_opt.alpha2 = alpha2; extra_opt.theta = theta;
        if mod(L,2)==0
            
            PSD = [fun(posf(1)),fun(posf(2:end-1))/2,...
                fun(posf(end)),fun(fliplr(posf(2:end-1)))/2]; %same way as components would be stored by using fft
            if debug_flag
                abs_freqs=[posf,fliplr(posf(2:end-1))];
                PSDpos = PSD(1:length(posf));
                PSDpos(2:end-1) = 2*PSDpos(2:end-1);
            end
            
        else
            
            PSD = [fun(posf(1)),fun(posf(2:end))/2,...
                fun(fliplr(posf(2:end)))/2];
            if debug_flag
                PSDpos = PSD(1:length(posf));
                PSDpos(2:end) = 2*PSDpos(2:end);
                abs_freqs=[posf,fliplr(posf(2:end))];
            end
        end
        
        s0dft = sqrt(L*PSD); %PSD to magnitude of rotor
        
        Phifk = 2*pi*rand(N,L);   %random phases
        Zfk = s0dft.*exp(complex(0,Phifk)); %add random phases while same frequency content maintained
        Iarray = ifft(Zfk,[],2,'symmetric'); %symmetric to get real output however also some artefacts but not important here (see calcgetfourier signal)
        Iarray = Iarray./max(abs(Iarray),[],2);
        
        if debug_flag
            %create plots in case debug needed
            for iN=1:N
                
                [N_p4(iN,:),N_f4(iN,:)] = genPSD(Iarray(iN,:),Tsim,Fs,'db');
            end
            iNs = ranperm(N,min(10,N));
            figure
            subplot(2,2,[1,2])
            for iN=iNs
                plot(Tsim,Iarray(iN,:),'displayName','Iarray_temp')
                hold on
            end
            legend('show')
            ylim([-1.5,1.5])
            subplot(2,2,3)
            plot(posf,10*log10(PSDpos),'displayName','fft:Mothersignal')
            
            set(gca,'xscale','log')
            legend('show')
            subplot(2,2,4)
            for iN=iNs
                plot(N_f4(iN,:),N_p4(iN,:),'displayName','Mothersignal+random phases')
                hold on
            end
            hold off
            
        end
        
    case lower('artificialPSD_BC')
        % the currents are derived from a aritificial PSD  based on
        % observed trends in literature for low frequency(<100Hz) and model
        % ouputs (Neuron) High frequency. A best case sprectrum is created, where a third breakpoint is added in order to drive frequency content to zero at highes >1e4 f.
        % This way minimizes interpolation errors. see An_spec script for
        % more details
        L = length(Tsim);
        Fs = 1/(Tsim(2)-Tsim(1));
        dF = Fs/L;
        posf = 0:dF:Fs/2; 
        
        wn = 14;
        theta = 0.7*rand()+0.3;         %max is reached between wn and theta*wn
        w2 = 10^3;                      %breakpoint higher powerlaw
        alpha1 = rand()+2;
        alpha2 = 2*rand()+2;
        alpha3 = 10;
        w3 = 10^4;
        %noise2 =@(x) 10.^(randn(1,length(x)));  %different noise methods
        %were tested see script An_spec
        noise3 = @(x) lognrnd(0,0.5,1,length(x));
        y = @(x) (1+x/(theta*wn))./((1+(x/wn).^alpha1).*(1+(x/w2).^alpha2).*(1+(x/w3).^alpha3));
        fun = @(x) y(x).*noise3(x); %this function determines artificial PSD right hand sided
        
        %store parameters
        extra_opt.alpha1 = alpha1; extra_opt.alpha2 = alpha2; extra_opt.theta = theta;
        if mod(L,2)==0
            
            PSD = [fun(posf(1)),fun(posf(2:end-1))/2,...
                fun(posf(end)),fun(fliplr(posf(2:end-1)))/2];
            if debug_flag
                abs_freqs=[posf,fliplr(posf(2:end-1))];
                PSDpos = PSD(1:length(posf));
                PSDpos(2:end-1) = 2*PSDpos(2:end-1);
            end
            
        else
            
            PSD = [fun(posf(1)),fun(posf(2:end))/2,...
                fun(fliplr(posf(2:end)))/2];
            if debug_flag
                PSDpos = PSD(1:length(posf));
                PSDpos(2:end) = 2*PSDpos(2:end);
                abs_freqs=[posf,fliplr(posf(2:end))];
            end
        end
        
        s0dft = sqrt(L*PSD); %PSD to magnitude of rotor
        
        Phifk = 2*pi*rand(N,L);   %random phases
        Zfk = s0dft.*exp(complex(0,Phifk)); %add random phases while same frequency content maintained
        Iarray = ifft(Zfk,[],2,'symmetric'); %symmetric to get real output however also some artefacts but not important here (see calcgetfourier signal)
        Iarray = Iarray./max(abs(Iarray),[],2);
        
        if debug_flag
            for iN=1:N
                
                [N_p4(iN,:),N_f4(iN,:)] = genPSD(Iarray(iN,:),Tsim,Fs,'db');
            end
            iNs = randi(N,1,min(10,N));
            figure
            subplot(2,2,[1,2])
            for iN=iNs
                plot(Tsim,Iarray(iN,:),'displayName','Iarray_temp')
                hold on
            end
            legend('show')
            ylim([-1.5,1.5])
            subplot(2,2,3)
            plot(posf,10*log10(PSDpos),'displayName','fft:Mothersignal')
            
            set(gca,'xscale','log')
            legend('show')
            subplot(2,2,4)
            for iN=iNs
                plot(N_f4(iN,:),N_p4(iN,:),'displayName','Mothersignal+random phases')
                hold on
            end
            hold off
            
        end
    case lower('artificialPSD_BC_val')
        % validation diffeernt methods ==> randomness excluded
        % the currents are derived from a aritificial PSD  based on
        % observed trends in literature for low frequency(<100Hz) and model
        % ouputs (Neuron) High frequency. A best case sprectrum is created, where a third breakpoint is added in order to drive frequency content to zero at highes >1e4 f.
        % This way minimizes interpolation errors. see An_spec script for
        % more details
        
        L = length(Tsim);
        Fs = 1/(Tsim(2)-Tsim(1));
        dF = Fs/L;
        posf = 0:dF:Fs/2;
        
        wn = 14;
        theta = 0.7*0.5+0.3;         %max is reached between wn and theta*wn
        w2 = 10^3;                      %breakpoint higher powerlaw
        alpha1 = 0.5+2;
        alpha2 = 2*0.5+2;
        alpha3 = 10;
        w3 = 10^4;
        %noise2 =@(x) 10.^(randn(1,length(x)));  %different noise methods
        %were tested see script An_spec
        %noise3 = @(x) lognrnd(0,0.5,1,length(x));
        fun = @(x) (1+x/(theta*wn))./((1+(x/wn).^alpha1).*(1+(x/w2).^alpha2).*(1+(x/w3).^alpha3));
        %fun = @(x) y(x).*noise3(x); %this function determines artificial PSD right hand sided
        
        
        %store parameters
        extra_opt.alpha1 = alpha1; extra_opt.alpha2 = alpha2; extra_opt.theta = theta;
        if mod(L,2)==0
            
            PSD = [fun(posf(1)),fun(posf(2:end-1))/2,...
                fun(posf(end)),fun(fliplr(posf(2:end-1)))/2];
            if debug_flag
                abs_freqs=[posf,fliplr(posf(2:end-1))];
                PSDpos = PSD(1:length(posf));
                PSDpos(2:end-1) = 2*PSDpos(2:end-1);
            end
            
        else
            
            PSD = [fun(posf(1)),fun(posf(2:end))/2,...
                fun(fliplr(posf(2:end)))/2];
            if debug_flag
                PSDpos = PSD(1:length(posf));
                PSDpos(2:end) = 2*PSDpos(2:end);
                abs_freqs=[posf,fliplr(posf(2:end))];
            end
        end
        
        s0dft = sqrt(L*PSD); %PSD to magnitude of rotor
        M = ([linspace(1,2,N)]'.*[linspace(1,2,L)]-ones(N,L))./3;%rand(N,L)
        Phifk = 2*pi*M;   %random phases
        Zfk = s0dft.*exp(complex(0,Phifk)); %add random phases while same frequency content maintained
        Iarray = ifft(Zfk,[],2,'symmetric'); %symmetric to get real output however also some artefacts but not important here (see calcgetfourier signal)
        Iarray = Iarray./max(abs(Iarray),[],2);
        
        if debug_flag
            for iN=1:N
                
                [N_p4(iN,:),N_f4(iN,:)] = genPSD(Iarray(iN,:),Tsim,Fs,'db');
            end
            iNs = randi(N,1,min(10,N));
            figure
            subplot(2,2,[1,2])
            for iN=iNs
                plot(Tsim,Iarray(iN,:),'displayName','Iarray_temp')
                hold on
            end
            legend('show')
            ylim([-1.5,1.5])
            subplot(2,2,3)
            plot(posf,10*log10(PSDpos),'displayName','fft:Mothersignal')
            
            set(gca,'xscale','log')
            legend('show')
            subplot(2,2,4)
            for iN=iNs
                plot(N_f4(iN,:),N_p4(iN,:),'displayName','Mothersignal+random phases')
                hold on
            end
            hold off
            
        end
    case lower('artificialPSD_WC')
        % the currents are derived from a aritificial PSD  based on
        % observed trends in literature for low frequency(<100Hz) and model
        % ouputs (Neuron) High frequency. A worst case sprectrum is
        % created. here the low frequency powerlw is extrapolated
        % throughout whole frequency domain
        % This way minimizes interpolation errors. see An_spec script for
        % more details
        L = length(Tsim);
        Fs = 1/(Tsim(2)-Tsim(1));
        dF = Fs/L;
        posf = 0:dF:Fs/2; 
        
        wn = 14;
        theta = 0.7*rand()+0.3;         %max is reached between wn and theta*wn
        alpha1 = rand()+2;

        %noise2 =@(x) 10.^(randn(1,length(x)));  %different noise methods
        %were tested see script An_spec
        noise3 = @(x) lognrnd(0,0.5,1,length(x));
        y = @(x) (1+x/(theta*wn))./((1+(x/wn).^alpha1));
        fun = @(x) y(x).*noise3(x); %this function determines artificial PSD right hand sided
        %store parameters
        extra_opt.alpha1 = alpha1; extra_opt.theta = theta;
        if mod(L,2)==0
            
            PSD = [fun(posf(1)),fun(posf(2:end-1))/2,...
                fun(posf(end)),fun(fliplr(posf(2:end-1)))/2];
            if debug_flag
                abs_freqs=[posf,fliplr(posf(2:end-1))];
                PSDpos = PSD(1:length(posf));
                PSDpos(2:end-1) = 2*PSDpos(2:end-1);
            end
            
        else
            
            PSD = [fun(posf(1)),fun(posf(2:end))/2,...
                fun(fliplr(posf(2:end)))/2];
            if debug_flag
                PSDpos = PSD(1:length(posf));
                PSDpos(2:end) = 2*PSDpos(2:end);
                abs_freqs=[posf,fliplr(posf(2:end))];
            end
        end
        
        s0dft = sqrt(L*PSD); %PSD to magnitude of rotor
        
        Phifk = 2*pi*rand(N,L);   %random phases
        Zfk = s0dft.*exp(complex(0,Phifk)); %add random phases while same frequency content maintained
        Iarray = ifft(Zfk,[],2,'symmetric'); %symmetric to get real output however also some artefacts but not important here (see calcgetfourier signal)
        Iarray = Iarray./max(abs(Iarray),[],2);
        
        if debug_flag
            for iN=1:N
                
                [N_p4(iN,:),N_f4(iN,:)] = genPSD(Iarray(iN,:),Tsim,Fs,'db');
            end
            iNs = randi(N,1,min(10,N));
            figure
            subplot(2,2,[1,2])
            for iN=iNs
                plot(Tsim,Iarray(iN,:),'displayName','Iarray_temp')
                hold on
            end
            legend('show')
            ylim([-1.5,1.5])
            subplot(2,2,3)
            plot(posf,10*log10(PSDpos),'displayName','fft:Mothersignal')
            
            set(gca,'xscale','log')
            legend('show')
            subplot(2,2,4)
            for iN=iNs
                plot(N_f4(iN,:),N_p4(iN,:),'displayName','Mothersignal+random phases')
                hold on
            end
            hold off
            
        end
        
    case 'randomdelay'
        % a predefined signal (see Itime) with a random time delay is given
        % to each dipole
        rnumbers = randi(length(Tsim),1,N);
        dt = Tsim(2)-Tsim(1);
        if strcmpi(Itime_type,'old')
            Tend_temp = 0.029;
        elseif strcmpi(Itime_type,'hipPyr25us1')
            Tend_temp = 0.1;
        end
        %Tend_temp = ceil(Tsim(end)/0.03)*0.03; %multiple of 30ms why necessary?
        Tsim_temp = 0:dt:Tend_temp;
        Isingle = Itime(Tsim_temp,d,Itime_type);
        Isingle = Isingle/max(abs(Isingle));
        for i = N:-1:1
            I_temp = circshift(Isingle,rnumbers(i));
            Iarray(i,:) = I_temp(1:length(Tsim));
            if Display
                Pval = round((N-i)/N*100,0);
                if Pval >= Pflag
                    disp(['Calculating I(t): ',num2str(Pval),'%'])
                    Pflag = Pflag+10;
                end
            end
        end
    case 'randomi'
        % load from a database currents with similar psd but random phases
        % (slow, database too large)
        if Tsim(end)>0.029
            if Tsim>0.058
                error('to long input')
            else
                load_temp = load('100randomIlength58ms.mat');
                Imat = load_temp.signals;
                t_Imat = load_temp.t2;
            end
        else
            load_temp = load('200randomIlength29ms.mat');
            Imat = load_temp.signals;
            t_Imat = load_temp.t;
        end
        clear('load_temp')
        if isempty(rselect)
            rselect = randi(size(Imat,1),1,N);
        elseif length(rselect)~=N
            error('false input rselect')
        end
        for i=N:-1:1
            I_sel = Imat(rselect(i),:)./max(abs(Imat(rselect(i),:)));
            I_temp = interp1(t_Imat,I_sel,Tsim);
            Iarray(i,:) = I_temp;
            if Display
                Pval = round((N-i+1)/N*100,0);
                if Pval >= Pflag
                    disp(['Calculating I(t): ',num2str(Pval),'%'])
                    Pflag = Pflag+10;
                end
            end
        end
    case 'randominewgen'
        % created currents with similar psd from mother current Itime and
        % add random phases. (more interpolation errors than artificial
        % PSD)
        
        if strcmpi(Itime_type,'old')
            Tend = 0.029;
        elseif strcmpi(Itime_type,'hipPyr25us1')
            Tend = 0.1;
        end
        % init signal
        nTend = ceil(Tsim(end)/Tend);
        Tend = nTend*Tend;
        Fs = 2^(17+nextpow2(nTend))/Tend; %this way Fs is at least 4MHz
        t = 0:1/Fs:Tend;
        s0 = Itime(t,d,Itime_type);  %upsampled signal with linear interp. create high freaquency components are questionqble 1/f relationship
        % calculate PSD
        L = length(s0);
        s0dft = fft(s0);
        filterdt = 100*1e-6;
        windowSize = filterdt*Fs;
        b = (1/windowSize)*ones(1,int16(windowSize));
        if ~PC_flag
            for i=N:-1:1
                Phifk = 2*pi*rand(1,L);   %random phases
                Zfk = s0dft.*exp(complex(0,Phifk)); %add random phases while same frequency content maintained
                Iarray_temp = ifft(Zfk,'symmetric'); %symmetric to get real output however also some artefacts but not important here (see calcgetfourier signal)
                Iarray_temp = Iarray_temp./max(abs(Iarray_temp));
                Iarray(i,:) = interp1(t,Iarray_temp,Tsim);
                Iarray(i,:) = filter(b,1,Iarray(i,:));
                Iarray(i,:) = Iarray(i,:)./max(Iarray(i,:));
                if Display
                    Pval = round((N-i+1)/N*100,0);
                    if Pval >= Pflag
                        disp(['Calculating I(t): ',num2str(Pval),'%'])
                        Pflag = Pflag+10;
                    end
                end
            end
        else
            parfor i=1:N
                Phifk = 2*pi*rand(1,L);
                Zfk = s0dft.*exp(complex(0,Phifk));
                Iarray_temp = ifft(Zfk,'symmetric');
                Iarray_temp = Iarray_temp./max(abs(Iarray_temp));
                Iarray(i,:) = interp1(t,Iarray_temp,Tsim);
                Iarray(i,:) = filter(b,1,Iarray(i,:));
                Iarray(i,:) = Iarray(i,:)./max(Iarray(i,:));
            end
        end
        
        
        
    case 'randomirandomdelay'
        % load from a database currents with similar psd but random phases
        % (slow, database too large). Random time delay added
        if Tsim(end)>0.029
            if Tsim>0.058
                error('to long input')
            else
                load_temp = load('100randomIlength58ms.mat');
                Imat = load_temp.signals;
                t_Imat = load_temp.t2;
            end
        else
            load_temp = load('200randomIlength29ms.mat');
            Imat = load_temp.signals;
            t_Imat = load_temp.t;
        end
        clear('load_temp')
        if isempty(rselect)
            rselect = randi(size(Imat,1),1,N);
        elseif length(rselect)~=N
            error('false input rselect')
        end
        rnumbers = randi(length(Tsim),1,N);
        for i=N:-1:1
            I_sel = Imat(rselect(i),:)./max(abs(Imat(rselect(i),:)));
            I_temp = interp1(t_Imat,circshift(I_sel,rnumbers(i)),Tsim);
            Iarray(i,:) = I_temp;
            if Display
                Pval = round((N-i)/N*100,0);
                if Pval >= Pflag
                    disp(['Calculating I(t): ',num2str(Pval),'%'])
                    Pflag = Pflag+10;
                end
            end
        end
    
        
        
        
    case 'sync'
        %all dipoles have same current assigned
        Isingle = Itime(Tsim,d,Itime_type);
        Iarray = repmat(Isingle,N,1);
    otherwise
        error('wrong type input')
end

% the maximum amplitude follows certain distribution.
switch lower(Itype)
    case 'same'
        Iarray = I*Iarray;
    case 'val'
        Iarray = [linspace(0.5,2,size(Iarray,1))]'.*Iarray;
    case 'normal'
        Iarray = normrnd(I,sigmaI,size(Iarray,1),1).*Iarray;
    case 'lognormal'
        % I is mean of randn generated values, sigmaI is standard deviation
        % of lognormal distribution containing 95% in 2 decades ==> sigma =
        % 2*log(10)
        mu = log(I^2/sqrt(sigmaI^2+I^2));
        Iarray = lognrnd(mu,sigmaI,size(Iarray,1),1).*Iarray;
        Iarray = min(Iarray,1e3);
        
end
end