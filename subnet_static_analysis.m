clc;clear;close all
%% Static analysis: This function calculates the calibrated Hin, Hse and HB from fMRI time series 
N=400;
basedir="~/Programs/Hierarchical-module-analysis/data/signal/"
outputdir="~/Programs/Hierarchical-module-analysis/output/"

network_cuts = [30, 70, 93, 118, 131, 161, 200];
network_names = {"Vis", "SomMot", "DorsAttn", "SalVentAttn", "Limbic", "Cont", "Default"};

fmridir = strcat(basedir, sprintf("par%d/", N))

subjs = load(strcat(outputdir, sprintf("subjects_%d.mat", N))).subjs
% Remove empty entries (subjects that failed in static analysis)
subjs = subjs(~cellfun('isempty', subjs))

fmrits={};IN={};IM={};

network_cuts = [network_cuts, network_cuts + N / 2]
network_names = [
    cellfun(@(v) strcat("LH_", v), network_names, "UniformOutput", false),
    cellfun(@(v) strcat("RH_", v), network_names, "UniformOutput", false)
]

%mypool=parpool('local',24,'IdleTimeout',240);
for sub=1:length(subjs)
    subj = subjs{sub};
    subfile = strcat(fmridir, subj)
    subname = regexprep(subj, "^(sub\d+).*$", "$1");

    fmri = load(subfile).ROI_ts;
    parfor network_id=1:length(network_names)
        network_start = [];
        network_end = network_cuts(network_id);
        if network_id == 1
            network_start = 1;
        else
            network_start = network_cuts(network_id - 1);
        end

        subnet = fmri(network_start:network_end, network_start:network_end);
        fmrits{sub, network_id} = subnet';
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

        IN{sub, network_id} = Hin;
        IM{sub, network_id} = Hse;
    end
end
IN = cell2mat(IN);
IM = cell2mat(IM);

if isfolder(outputdir) == 0
    mkdir(outputdir);
end
save(strcat(outputdir, sprintf("origin_Hb_subnet_static_%d.mat", N)), 'IN', 'IM');

%% stable FC matrix for long-enough fMRI length
%sFC=corr(Long_fmri);
%% Calibrating the individual static segregation and integration component

Long_fmri = {}; sFC = {};
Hins = {}; Hses = {};
for i=1:7
    Long_fmri{i} = cell2mat(fmrits(:, 1));
    sFC{i} = corr(Long_fmri{i});
    [Hin, Hse] = Stable_correct(sFC{i}, IN(:, i), IM(:, i), length(sFC{i}));
    Hins{i} = Hin; Hses{i} = Hse;
end
save(strcat(outputdir, sprintf("corrected_Hb_subnet_static_%d.mat", N)), 'Hins', 'Hses');
