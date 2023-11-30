% print correct intensity and pressure fields
function printCalcveloc(Input,targetA)

close all
disp('start function')
if ischar(Input)    
Input_temp = load(['./Inputs/',Input,'.mat']);
Input = Input_temp.Input;
end
  
S4l_flag = Input.S4l_flag; S4l_fn = Input.Filename_S4l; S4l_loc = Input.Location_S4l;
Model = Input.Model;
ROI = Input.ROI;

switch lower(Model)
    case 'human_scalp'
        SolutionType = '4SphereS8.7R25';
        RPOI = 0.082;
        RBrain = 0.07;
        dRBrain = 0.005;
    case 'human_cortical'
        SolutionType = '4SphereS8.7R25';
        RPOI = 0.07;
        RBrain = 0.07;
        dRBrain = 0.005;
    case 'mouse_fair0'
        SolutionType = 'Mouse4Sphere~fair0';
        RPOI = 0.0059;
        RBrain = 0.0046;
        dRBrain = 0.0005;
    otherwise 
        error('false input model')        
end
PLOT_s4l = 1;




% load pressure field
S4l_input = struct();
[~, ~, ~, fus, maxv,~,~] = ...
    calcvelocity(S4l_loc,S4l_fn,RBrain-dRBrain,0,'snapshot',0,'rotate',0','totalperiods',50,'inittime',48);
S4l_input.loc.unit = 'mm';    
vib_Amp = maxv/(fus*2*pi);
Scale = targetA/vib_Amp;
[~, ~, ~, ~, maxv2,maxv_loc,~] = ...
    calcvelocity(S4l_loc,S4l_fn,RBrain-dRBrain,1,'scaleFactor',Scale,'snapshot',0,'rotate',0','totalperiods',50,'inittime',48);
new_vib_Amp = maxv2/(2*pi*fus);
if abs(targetA-new_vib_Amp)<targetA/100
    fprintf('good job')
else
    fprintf(':(')
end
pos_string = sprintf('%3.2f-',maxv_loc);
prefix = [Model,'_',ROI,'_',pos_string(1:end-1),'_',sprintf('%3.2e',targetA)];
suffix = {'pvs','pvdsl','pdI'};
for i=1:3
    figname = ['./Figs pressure fields/',prefix,'_',suffix{i},datestr(now,'ddmmyy'),'.fig'];
    f = figure(i);
    savefig(f,figname)
end
end