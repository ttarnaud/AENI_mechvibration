% Examples how to compile scripts (also applicable on HPC)

rmpath(genpath('D:\Users\rschoeters\Documents\MATLAB\Interestingfunctions'));
addpath('D:\Users\rschoeters\Documents\MATLAB\Interestingfunctions');

mcc -mv calcStrengthCurve_mosc_DBvsnormal.m -I ./DataBases -I ./Functions -I ./Inputs -I ./PressureFields...
    -I ./Functions/Subfunctions -I ./Functions/Dipole -d ./Compilations

%%
rmpath(genpath('D:\Users\rschoeters\Documents\MATLAB\Interestingfunctions'));
addpath('D:\Users\rschoeters\Documents\MATLAB\Interestingfunctions');

mcc -mv srun_cSC_mosc_DBvsNormal.m -I ./DataBases -I ./Functions -I ./Inputs -I ./PressureFields...
    -I ./Functions/Subfunctions -I ./Functions/Dipole -d ./Compilations
%%
rmpath(genpath('D:\Users\rschoeters\Documents\MATLAB\Interestingfunctions'));
addpath('D:\Users\rschoeters\Documents\MATLAB\Interestingfunctions');

mcc -mv run_cSC_mosc_DBvsNormal.m -I ./DataBases -I ./Functions -I ./Inputs -I ./PressureFields...
    -I ./Functions/Subfunctions -I ./Functions/Dipole -d ./Compilations

%%
rmpath(genpath('D:\Users\rschoeters\Documents\MATLAB\Interestingfunctions'));
addpath('D:\Users\rschoeters\Documents\MATLAB\Interestingfunctions');

mcc -mv run_cSC_mosc_DBfvsNormal.m -I ./DataBases -I ./Functions -I ./Inputs -I ./PressureFields...
    -I ./Functions/Subfunctions -I ./Functions/Dipole -d ./Compilations
%%
rmpath(genpath('D:\Users\rschoeters\Documents\MATLAB\Interestingfunctions'));
addpath('D:\Users\rschoeters\Documents\MATLAB\Interestingfunctions');

mcc -mv srun_cSC_mosc_DBfvsNormal.m -I ./DataBases -I ./Functions -I ./Inputs -I ./PressureFields...
    -I ./Functions/Subfunctions -I ./Functions/Dipole -d ./Compilations

%%
% compile script

rmpath(genpath('D:\Users\rschoeters\Documents\MATLAB\Interestingfunctions'));
addpath('D:\Users\rschoeters\Documents\MATLAB\Interestingfunctions');

mcc -mv ./Functions/calcErrorGrid.m -I ./Functions -I ./Inputs -I ./Functions/Subfunctions ...
    -I ./Functions/Dipole -d ./Compilations