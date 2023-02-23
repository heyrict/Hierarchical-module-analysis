clc;clear;close all
%% Static analysis: This function calculates the calibrated Hin, Hse and HB from fMRI time series 
N=400;
basedir="~/Programs/Hierarchical-module-analysis/data/signal/"
outputdir="~/Programs/Hierarchical-module-analysis/output/"

fmridir = strcat(basedir, sprintf("par%d/", N))

subjs = load(strcat(outputdir, sprintf("subjects_%d.mat", N))).subjs
% Remove empty entries (subjects that failed in static analysis)
subjs = subjs(~cellfun('isempty', subjs))

fmrits={};IN={};IM={};

%mypool=parpool('local',24,'IdleTimeout',240);
parfor sub=1:length(subjs)
    subj = subjs{sub};
    subfile = strcat(fmridir, subj)
    subname = regexprep(subj, "^(sub\d+).*$", "$1");

    fmri = load(subfile).ROI_ts;
    fmrits{sub} = fmri';
    %% individual static FC matrix and its hierarchical module partition
    FC=corr(fmri);
    Clus_num = [];
    Clus_size = [];
    mFC = []

    % Calculating H_i
    [Clus_num,Clus_size,mFC] = Functional_HP(FC,N);
    FC=(FC+FC')/2;
    FC(FC<0)=0;
    [FEC FE]=eig(FC);
    FE(FE<0)=0;
    FE=FE^2;%% using the squared Lambda
    p=zeros(1,N);
    for i=1:length(find(Clus_num<1))
          p(i)=sum(abs(Clus_size{i}-1/Clus_num(i)))/N;%% modular size correction
    end
    HF=fliplr(diag(FE)').*Clus_num.*(1-p);

    % Calculating parcel-based integration and segregation
    Hin = HF(1) .* FEC(1, :) .^ 2;
    Hse = zeros(N);
    for j=1:N
        Hse(j) = sum(HF .* FEC(j, :) .^ 2);
    end

    IN{sub} = Hin;
    IM{sub} = Hse;
end
Long_fmri = cell2mat(fmrits)';
IN = cell2mat(IN);
IM = cell2mat(IM);

if isfolder(outputdir) == 0
    mkdir(outputdir);
end
save(strcat(outputdir, sprintf("subjects_%d.mat", N)), 'subjs');
save(strcat(outputdir, sprintf("origin_Hb_parcel_static_%d.mat", N)), 'IN', 'IM');

%% stable FC matrix for long-enough fMRI length
sFC=corr(Long_fmri);
%% Calibrating the individual static segregation and integration component

Hins = {}; Hses = {};
[Hins, Hses] = Stable_correct(sFC, IN, IM, N);
save(strcat(outputdir, sprintf("corrected_Hb_parcel_static_%d.mat", N)), 'Hins', 'Hses');
