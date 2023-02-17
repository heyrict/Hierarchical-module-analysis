clc;clear;close all
%% This function generates the temporally dynamic Hin and Hse, and then calibrates them within individuals. 
N = 400;
SeqLen=220;width=32;step=1;TR=2
basedir="~/Programs/Hierarchical-module-analysis/data/signal/"
outputdir="~/Programs/Hierarchical-module-analysis/output/"

subjs = load(strcat(outputdir, sprintf("subjects_%d.mat", N))).subjs
% Remove empty entries (subjects that failed in static analysis)
subjs = subjs(~cellfun('isempty', subjs))

network_cuts = [30, 70, 93, 118, 131, 161, 200];
network_names = {"Vis", "SomMot", "DorsAttn", "SalVentAttn", "Limbic", "Cont", "Default"};

network_cuts = [network_cuts, network_cuts + N / 2]
network_names = [
    cellfun(@(v) strcat("LH_", v), network_names, "UniformOutput", false),
    cellfun(@(v) strcat("RH_", v), network_names, "UniformOutput", false)
]

fmridir = strcat(basedir, sprintf("par%d/", N))
N_sub=length(subjs);
%%%===============================================
% mypool=parpool('local',24,'IdleTimeout',240);
S=load(strcat(outputdir, sprintf('corrected_Hb_static_%d.mat', N)));
Z=[];Hins=[];Hses=[];cHins=[];cHses=[];

if isfolder(outputdir) == 0
    mkdir(outputdir);
end

for sub=1:N_sub;
    % Retrieve subject id from string
    subj = subjs{sub};
    subfile = strcat(fmridir, subj)
    subname = regexprep(subj, "^(sub\d+).*$", "$1");

    fmri = load(subfile).ROI_ts;
    IN={};IM={};
    for network_id=1:length(network_names)
        network_start = [];
        network_end = network_cuts(network_id);
        if network_id == 1
            network_start = 1;
        else
            network_start = network_cuts(network_id - 1);
        end

        subnet = fmri(:, network_start:network_end);
        num_loops = floor((SeqLen - width) / step) + 1;
        parfor l=1:num_loops
            t = (l - 1) * step + 1;
            subdata=subnet(t:t+width-1, :);
            FC=corr(subdata);
            try
                [Clus_num,Clus_size] = Functional_HP(FC, network_end - network_start + 1);
            catch
                warning(strcat("Error executing Functional_HP, Filename: ", subfile))
                continue
            end
            [Hin,Hse] =Balance(FC,network_end - network_start + 1,Clus_size,Clus_num);
            IN{l, network_id} = Hin;
            IM{l, network_id} = Hse;
        end
    end
    IN = cell2mat(IN);
    IM = cell2mat(IM);
    Hins=[Hins;IN];
    Hses=[Hses;IM];
    [Hin] = individual_correction(IN,S.Hin(sub));
    [Hse] = individual_correction(IM,S.Hse(sub));
    cHins=[cHins;Hin];
    cHses=[cHses;Hse];
    [Fre,DIn,DSe,In_time,Se_time] = Flexible(Hin-Hse, TR)%% calculating dynamic measures
    Z=[Z;Fre,DIn,DSe,In_time,Se_time];
end

save(strcat(outputdir, sprintf("Hb_dynamic_%d_w%d_s%d.mat", N, width, step)), 'Hins', 'Hses', 'cHins', 'cHses');
save(strcat(outputdir, sprintf("measures_dynamic_%d_w%d_s%d.mat", N, width, step)), 'Z', 'SeqLen', 'width', 'step', 'TR');
