clc;clear;close all
%% Static analysis: This function calculates the calibrated Hin, Hse and HB from fMRI time series 
N=100;
basedir="~/Programs/Hierarchical-module-analysis/data/signal/"
outputdir="~/Programs/Hierarchical-module-analysis/output/"

fmridir = strcat(basedir, sprintf("par%d/", N))
subjs = dir(fmridir);
N_sub=length(subjs);
truesubjs = {};
for sub=1:N_sub % Loop through all files in directory
    % Retrieve subject id from string
    subj = subjs(sub).name;
    if !length(regexp(subj, "^sub.*\\.mat$", "once"))
        continue;
    end
    truesubjs = [truesubjs, subj];
end
subjs = truesubjs;

fmrits={};IN={};IM={};subjs_processed = {};

%mypool=parpool('local',24,'IdleTimeout',240);
parfor sub=1:length(subjs)
    subj = subjs{sub};
    subfile = strcat(fmridir, subj)
    subname = regexprep(subj, "^(sub\\d+).*$", "$1");

    fmri = load(subfile).ROI_ts;
    %% individual static FC matrix and its hierarchical module partition
    FC=corr(fmri);
	Clus_num = [];
	Clus_size = [];
	mFC = []
    try
		[Clus_num,Clus_size,mFC] = Functional_HP(FC,N);
	catch
		warning(strcat("Error executing Functional_HP, Filename: ", subfile))
		continue
	end
    %parsave(subname,Clus_size,Clus_num,'_Clus.mat')
    %parsave(subname,FC,mFC,'_FC.mat')
    %% individual static integration component Hin and segragtion component Hse
    [Hin,Hse] =Balance(FC,N,Clus_size,Clus_num);

	subjs_processed{sub} = subj
    fmrits{sub} = fmri';
    IN{sub} = Hin;
    IM{sub} = Hse;
end
Long_fmri = cell2mat(fmrits)';
IN = cell2mat(IN);
IM = cell2mat(IM);

if isfolder(outputdir) == 0
    mkdir(outputdir);
end
subjs = subjs_processed;
save(strcat(outputdir, sprintf("subjects_%d.mat", N)), 'subjs');
save(strcat(outputdir, sprintf("origin_Hb_static_%d.mat", N)), 'IN', 'IM');

%% stable FC matrix for long-enough fMRI length
sFC=corr(Long_fmri);
%% Calibrating the individual static segregation and integration component
[Hin,Hse] = Stable_correct(sFC,IN,IM,N);
save(strcat(outputdir, sprintf("corrected_Hb_static_%d.mat", N)), 'Hin', 'Hse');
