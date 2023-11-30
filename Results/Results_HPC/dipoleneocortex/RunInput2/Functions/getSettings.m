function [Options,RSphere,RPOI,RTICKS,Rvals,Angles] = getSettings(ModelType,resolution)

fDependence = 0;
tScalp = 0.007;
tSkull = 0.005;
sigma = 0.33;
Val_flag = 0;
tAir = 0.005;
dAngle = pi/(10*resolution);
Angles = 0:dAngle:pi;
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
        RPOI = RSphere+0.017;
        SolutionType = '3Sphere';
        RatioSkT = 1/25;
        Rvals = linspace(0,RSphere-0.002,resolution+1);
        RTICKS = unique([0:0.01:Rvals(end),Rvals(end),fliplr(RPOI:-0.01:Rvals(end))],'stable');
        tScalp = 0.012;
    case lower('4SphereS8.7R~f_{VAL}')
        RSphere = 0.07;
        RPOI = RSphere+0.017;
        SolutionType = '4Sphere';
        RatioSkT = 1/25;
        Rvals = linspace(0,RSphere-0.002,resolution+1);
        RTICKS = unique([0:0.01:Rvals(end),Rvals(end),fliplr(RPOI:-0.01:Rvals(end))],'stable');
        fDependence = 1;
        sigma = [0.33,RatioSkT*0.33,0.33,0.33];
        Val_flag = 1;
    case lower('4SphereS8.7R25_{VAL}')
        RSphere = 0.07;
        RPOI = RSphere+0.017;
        SolutionType = '4Sphere';
        RatioSkT = 1/25;
        Rvals = linspace(0,RSphere-0.002,resolution+1);
        RTICKS = unique([0:0.01:Rvals(end),Rvals(end),fliplr(RPOI:-0.01:Rvals(end))],'stable');
        sigma = [0.33,RatioSkT*0.33,0.33,0.33];
    case lower('4SphereS8.7R25_{VAL}')
        RSphere = 0.07;
        RPOI = RSphere+0.017;
        SolutionType = '4Sphere';
        RatioSkT = 1/25;
        Rvals = linspace(0,RSphere-0.002,resolution+1);
        RTICKS = unique([0:0.01:Rvals(end),Rvals(end),fliplr(RPOI:-0.01:Rvals(end))],'stable');
        sigma = [0.33,RatioSkT*0.33,0.33,0.33];
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