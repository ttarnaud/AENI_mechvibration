function out = alphaTrainFun(t,varargin)
%Create input based on train of alpha's
% Required input is the time t
% - Either a string pointing to file that is loaded containing structure with
% 3 fields: Taus,AlphaDelay and Amps
% or three additional inputs:
%     Taus = varargin{1};
%    AlphaDelay = varargin{2};
%    Amps = varargin{3};

if nargin ==2
    filename = varargin{1};
    try 
        load(filename);
    catch
        load(['./Inputs/',filename]);
    end
        Taus = aTInput.Taus;
        AlphaDelay = aTInput.AlphaDelay;
        Amps = aTInput.Amps;
elseif nargin ==4
    Taus = varargin{1};
    AlphaDelay = varargin{2};
    Amps = varargin{3};
else
end
Alphafun = @(t,Tau,t0,A) A*double((t-t0)>=0).*((t-t0)./Tau).*exp(1-(t-t0)./Tau)+0; %alpha function
Ifun =@(t) sum(cell2mat(arrayfun(@(Tau,AlphaDelay,Amp) Alphafun(t,Tau,AlphaDelay,Amp),Taus,AlphaDelay,Amps,'UniformOutput',0)'),1);
I = Ifun(t);
out = I/(I(abs(I)==max(abs(I))));
end