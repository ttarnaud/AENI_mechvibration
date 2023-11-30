difference with version 2: 

%Inputs
intm = load('./Results/SCs_all_cEG_kAus_nos4l_v2.mat');
requiredAmpAll = load('D:\users\rschoeters\Documents\Imec USEEG\Matlab\Results\SCs_all_cEG_Amp_nos4l_fkApdpset_v2.mat', 'requiredAmpAll');
% narrow down Ampgrid
requiredAmpAll = requiredAmpAll.requiredAmpAll(:,idxdps);
min_rAA = floor(log10(min(requiredAmpAll,[],2)))-1;
max_rAA = ceil(log10(max(requiredAmpAll,[],2)))+1;

%% ----------difference is here below----------------
if strcmpi(dpI_time,'artificialPSD_fixed3')
    max_rAA = min_rAA;
    min_rAA = min_rAA-3;
elseif strcmpi(dpI_time,'artificialPSD_fixed2')
    max_rAA = max_rAA-1;
    min_rAA = min_rAA-1;
end
------------------------------------------------------

Amp = arrayfun(@(a,b)10.^(a:b),min_rAA,max_rAA,'UniformOutput',false);