clc;clear;close all
%% Static analysis: This function calculates the calibrated Hin, Hse and HB from fMRI time series 
N=200;
%basedir="/run/media/heyrict/wd22a/research/blackblood/signal/"
basedir="/tmp/data/signal/"
outputdir="/tmp/data/output/"

fmridir = strcat(basedir, "par", sprintf("%d", N), "/")
subjs = dir(fmridir);
N_sub=length(subjs);
truesubjs = {}
for sub=1:N_sub % Loop through all files in directory
    % Retrieve subject id from string
    subj = subjs(sub).name;
    if !length(regexp(subj, "^sub.*\\.mat$", "once"))
        continue;
    end
    truesubjs = [truesubjs, subj]
end
subjs = truesubjs

if !isfolder(outputdir)
    mkdir(outputdir);
end
save("-7", strcat(outputdir, "subjects_static.mat"), 'subjs');

fmrits={};IN={};IM={};

%mypool=parpool('local',24,'IdleTimeout',240);
parfor sub=1:length(subjs)
    subj = char(subjs(sub));
    subfile = strcat(fmridir, subj)
    subname = regexprep(subj, "^(sub\\d+).*$", "$1");

    fmri = load(subfile).ROI_ts;
    %% individual static FC matrix and its hierarchical module partition
    FC=corr(fmri);
    [Clus_num,Clus_size,mFC] = Functional_HP(FC,N);
    %parsave(subname,Clus_size,Clus_num,'_Clus.mat')
    %parsave(subname,FC,mFC,'_FC.mat')
    %% individual static integration component Hin and segragtion component Hse
    [Hin,Hse] =Balance(FC,N,Clus_size,Clus_num);

    fmrits(sub) = fmri';
    IN(sub) = Hin;
    IM(sub) = Hse;
end
Long_fmri = cell2mat(fmrits)';
IN = cell2mat(IN);
IM = cell2mat(IM);

save("-7", strcat(outputdir, "origin_Hb_static.mat"), 'IN', 'IM');

%% stable FC matrix for long-enough fMRI length
sFC=corr(Long_fmri);
%% Calibrating the individual static segregation and integration component
[Hin,Hse] = Stable_correct(sFC,IN,IM,N);
save("-7", strcat(outputdir, "corrected_Hb_static.mat"), 'Hin', 'Hse');
