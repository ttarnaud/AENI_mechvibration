function [Options,RSphere,RPOI,RTICKS,Rvals,Angles] = getSettings(ModelType,resolution)
%predefined geometric configurations of headmodels. Better to do this way
%in order to have consitency cross simulations. 
%Output contains a 
%       Options cell {'RatioSkT',RatioSkT,'RSphere',RSphere,...
%            'SolutionType',SolutionType,'fDependence',fDependence,'tSkull',tSkull,...
%           'tScalp',tScalp,'sigma',sigma,'Validation',Val_flag,'tAir',tAir}
%       Rsphere is size of brain in [m]
%       RPOI position of POIs = outer sphere
%       RTICKS: used in certain plots
%       Rvals: Rvals to go over when displacing DOI (used in first
%       simulations)
%       Angles: similar to RVals only the coelevetion here
%Input: Modeltype see below different cases
%       resolution: used to determine Rvals and Angles resolution

fDependence = 0;     %frequency dependend conductivity flag
tScalp = 0.007;      % standard thickness scalp in human head model [m]
tSkull = 0.005;      % standard thickness skull in human head model [m]
sigma = 0.33;        % conductivity of brain tissue
Val_flag = 0;        %Validation flag
tAir = 0.005;        % thickness of air layer
dAngle = pi/(10*resolution);    %steps for angles
Angles = 0:dAngle:pi;        

%solution type is important! defines which simulaiton model is used
switch lower(ModelType)
    case lower('ClosedBounded')
        RSphere = 0.07;
        RPOI = 0.07;
        RatioSkT = 1;
        SolutionType = 'ClosedBounded';
        Rvals = linspace(0,RPOI-0.002,resolution+1);
        RTICKS = unique([0:0.01:Rvals(end),Rvals(end),fliplr(RPOI:-0.01:Rvals(end))],'stable');
    case lower('ClosedBounded2')
        RSphere = 0.07;
        RPOI = 0.07;
        RatioSkT = 1;
        SolutionType = 'ClosedBounded';
        Rvals = linspace(0,0.058-0.002,resolution+1);
        RTICKS = unique([0:0.01:Rvals(end),Rvals(end),fliplr(RPOI:-0.01:Rvals(end))],'stable');
    case lower('4SphereasClosedBounded2')
        RSphere = 0.058;
        RPOI = 0.07;
        SolutionType = '4Sphere';
        RatioSkT = 1;
        Rvals = linspace(0,0.058-0.002,resolution+1);
        RTICKS = unique([0:0.01:Rvals(end),Rvals(end),fliplr(RPOI:-0.01:Rvals(end))],'stable');
        sigma = [0.33,RatioSkT*0.33,0.33,2*pi*8.85*10^-12*10^5;];
    case lower('3SphereS7R1Rval6.8')
        RSphere = 0.058;
        RPOI = RSphere+0.012;
        SolutionType = '3Sphere';
        RatioSkT = 1;
        Rvals = linspace(0,RPOI-0.002,resolution+1);
        RTICKS = unique([0:0.01:Rvals(end),Rvals(end),fliplr(RPOI:-0.01:Rvals(end))],'stable');
    case lower('3SphereS7R1')
        RSphere = 0.058;
        RPOI = RSphere+0.012;
        SolutionType = '3Sphere';
        RatioSkT = 1;
        Rvals = linspace(0,RSphere-0.002,resolution+1);
        RTICKS = unique([0:0.01:Rvals(end),Rvals(end),fliplr(RPOI:-0.01:Rvals(end))],'stable');
    case lower('3SphereS8.2R1')
        RSphere = 0.07;
        RPOI = RSphere+0.012;
        SolutionType = '3Sphere';
        RatioSkT = 1;
        Rvals = linspace(0,RSphere-0.002,resolution+1);
        RTICKS = unique([0:0.01:Rvals(end),Rvals(end),fliplr(RPOI:-0.01:Rvals(end))],'stable');
    case lower('3SphereS7R25')
        RSphere = 0.058;
        RPOI = RSphere+0.012;
        SolutionType = '3Sphere';
        RatioSkT = 1/25;
        Rvals = linspace(0,RSphere-0.002,resolution+1);
        RTICKS = unique([0:0.01:Rvals(end),Rvals(end),fliplr(RPOI:-0.01:Rvals(end))],'stable');
    case lower('3SphereS8.2R25')
        RSphere = 0.07;
        RPOI = RSphere+0.012;
        SolutionType = '3Sphere';
        RatioSkT = 1/25;
        Rvals = linspace(0,RSphere-0.002,resolution+1);
        RTICKS = unique([0:0.01:Rvals(end),Rvals(end),fliplr(RPOI:-0.01:Rvals(end))],'stable');
        tAir = 0;
    case lower('3SphereS8.2R~f')
        RSphere = 0.07;
        RPOI = RSphere+0.012;
        SolutionType = '3Sphere';
        RatioSkT = 1/25;
        Rvals = linspace(0,RSphere-0.002,resolution+1);
        RTICKS = unique([0:0.01:Rvals(end),Rvals(end),fliplr(RPOI:-0.01:Rvals(end))],'stable');
        fDependence = 1;
    case lower('4SphereS8.7R25')
        RSphere = 0.07;
        RPOI = RSphere+0.017;
        SolutionType = '4Sphere';
        RatioSkT = 1/25;
        Rvals = linspace(0,RSphere-0.002,resolution+1);
        RTICKS = unique([0:0.01:Rvals(end),Rvals(end),fliplr(RPOI:-0.01:Rvals(end))],'stable');
    case lower('4SphereS8.7R~f')
        tAir = 0.005;
        RSphere = 0.07;
        RPOI = RSphere+0.012+tAir;
        SolutionType = '4Sphere';
        RatioSkT = 1/25;
        Rvals = linspace(0,RSphere-0.002,resolution+1);
        RTICKS = unique([0:0.01:Rvals(end),Rvals(end),fliplr(RPOI:-0.01:Rvals(end))],'stable');
        fDependence = 1;
    case lower('3SphereS8.7R25_{VAL}')
        RSphere = 0.07;
        RPOI = RSphere+0.012;
        SolutionType = '3Sphere';
        RatioSkT = 1/25;
        Rvals = linspace(0,RSphere-0.002,resolution+1);
        RTICKS = unique([0:0.01:Rvals(end),Rvals(end),fliplr(RPOI:-0.01:Rvals(end))],'stable');
        tScalp = 0.007;
    case lower('4SphereS8.7R~f_{VAL}')
        RSphere = 0.07;
        RPOI = RSphere+0.012;
        SolutionType = '4Sphere';
        RatioSkT = 1/25;
        Rvals = linspace(0,RSphere-0.002,resolution+1);
        RTICKS = unique([0:0.01:Rvals(end),Rvals(end),fliplr(RPOI:-0.01:Rvals(end))],'stable');
        fDependence = 1;
        sigma = [0.33,RatioSkT*0.33,0.33,0];
        Val_flag = 1;
    case lower('4SphereS8.7R25_{VAL}')
        RSphere = 0.07;
        RPOI = RSphere+0.012;
        SolutionType = '4Sphere';
        RatioSkT = 1/25;
        Rvals = linspace(0,RSphere-0.002,resolution+1);
        RTICKS = unique([0:0.01:Rvals(end),Rvals(end),fliplr(RPOI:-0.01:Rvals(end))],'stable');
        sigma = [0.33,RatioSkT*0.33,0.33,0];
    case lower('Mouse4Sphere~fair0')
        tAir = 1e-20;
        RSphere = 0.0046;
        tSkull = 0.0006;
        tScalp = 0.0007;
        RPOI = RSphere+tSkull+tScalp+tAir;
        SolutionType = '4Sphere';
        RatioSkT = 1/25;
        Rvals = linspace(0,RSphere-0.0005,resolution+1);
        RTICKS = unique([0:0.01:Rvals(end),Rvals(end),fliplr(RPOI:-0.01:Rvals(end))],'stable');
    case lower('Mouse4Sphere~fair1')
        tAir = 0.001;
        RSphere = 0.0046;
        tSkull = 0.0006;
        tScalp = 0.0007;
        RPOI = RSphere+tSkull+tScalp+tAir;
        SolutionType = '4Sphere';
        RatioSkT = 1/25;
        Rvals = linspace(0,RSphere-0.0005,resolution+1);
        RTICKS = unique([0:0.01:Rvals(end),Rvals(end),fliplr(RPOI:-0.01:Rvals(end))],'stable');
    otherwise
        error('wrong input')
end

Options = {'RatioSkT',RatioSkT,'RSphere',RSphere,...
            'SolutionType',SolutionType,'fDependence',fDependence,'tSkull',tSkull,...
            'tScalp',tScalp,'sigma',sigma,'Validation',Val_flag,'tAir',tAir};

end