function Iarray = getIarray(type,Itype,Tsim,N,I,d,sigmaI,rselect)
Display = 0;
Pflag = 0;
Iarray = zeros(N,length(Tsim));
PC_flag = 1 && ~Display;

switch lower(type)
    case 'randomdelay'
        rnumbers = randi(length(Tsim),1,N);
        dt = Tsim(2)-Tsim(1);
        Tend_temp = ceil(Tsim(end)/0.03)*0.03;
        Tsim_temp = 0:dt:Tend_temp;
        Isingle = Itime(Tsim_temp,d);
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
        % init signal
        Tend=0.029;
        Fs = 2^17/Tend;
        t = 0:1/Fs:Tend-1/Fs;
        s0 = Itime(t,500e-6);
        % calculate PSD
        L = length(s0);
        s0dft = fft(s0);
        filterdt = 100*1e-6;
        windowSize = filterdt*Fs;
        b = (1/windowSize)*ones(1,int16(windowSize));
        if ~PC_flag
            for i=N:-1:1
                Phifk = 2*pi*rand(1,L);
                Zfk = s0dft.*exp(complex(0,Phifk));
                Iarray_temp = ifft(Zfk,'symmetric');
                Iarray_temp = Iarray_temp./max(abs(Iarray_temp));
                Iarray(i,:) = interp1(t,Iarray_temp,Tsim);
                Iarray(i,:) = filter(b,1,Iarray(i,:));
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
            end
        end
        
        
        
    case 'randomirandomdelay'
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
        Isingle = Itime(Tsim,d);
        Iarray = repmat(Isingle,N,1);
    otherwise
        error('wrong type input')
end
switch lower(Itype)
    case 'same'
        Iarray = I*Iarray;
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