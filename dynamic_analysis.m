clc;clear;close all
%% This function generates the temporally dynamic Hin and Hse, and then calibrates them within individuals. 
N = 200;
basedir="/tmp/data/signal/"
outputdir="/tmp/data/output/"

subjs = load(strcat(outputdir, "subjects_static.mat")).subjs

fmridir = strcat(basedir, "par", sprintf("%d", N), "/")
N_sub=length(subjs);
SeqLen=200;width=100;step=20;TR=2
subnames={};
%%%===============================================
% mypool=parpool('local',24,'IdleTimeout',240);
S=load(strcat(outputdir, 'corrected_Hb_static.mat'));
Z=[];Hins=[];Hses=[];cHins=[];cHses=[];

if !isfolder(outputdir)
    mkdir(outputdir);
end

for sub=1:N_sub;
    % Retrieve subject id from string
    subj = char(subjs(sub));
    if !length(regexp(subj, "^sub.*\\.mat$", "once"))
        continue;
    end
    subfile = strcat(fmridir, subj)
    subname = regexprep(subj, "^(sub\\d+).*$", "$1");
    subnames = [subnames, subname];

    fmri = load(subfile).ROI_ts;

    IN=[];IM=[];

    % Use for parallelized computing
    %function [Hin, Hse] = CalcHStep(t)
    %    subdata = fmri(t:t+width, :);
    %    FC = corr(subdata);
    %    [Clus_num, Clus_size] = Functional_HP(FC, N);
    %    [Hin, Hse] = Balance(FC, N, Clus_size, Clus_num);
    %end

    parfor t=1:step:SeqLen-width
        subdata=fmri(t:t+width, :);
        FC=corr(subdata);
        [Clus_num,Clus_size] = Functional_HP(FC,N);
        [Hin,Hse] =Balance(FC,N,Clus_size,Clus_num);
        IN=[IN,Hin];IM=[IM,Hse];
    end
    Hins=[Hins;IN];
    Hses=[Hses;IM];
    [Hin] = individual_correction(IN,S.Hin(sub));
    [Hse] = individual_correction(IM,S.Hse(sub));
    cHins=[cHins;Hin];
    cHses=[cHses;Hse];
    [Fre,DIn,DSe,In_time,Se_time] = Flexible(Hin-Hse, TR)%% calculating dynamic measures
    Z=[Z;Fre,DIn,DSe,In_time,Se_time];
end

%save("-7", strcat(outputdir, "subjects_dynamic.mat"), 'subnames');
save("-7", strcat(outputdir, "Hb_dynamic.mat"), 'Hins', 'Hses', 'cHins', 'cHses');
save("-7", strcat(outputdir, "measures_dynamic.mat"), 'Z', 'SeqLen', 'width', 'step', 'TR');
