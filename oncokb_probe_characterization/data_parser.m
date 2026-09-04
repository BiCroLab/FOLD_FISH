%A script to extract the features of FOLD-FISH probes targeting OncoKB and COSMIC genes
%designed by Raquel Silva for the FOLD-FISH manuscript
%
%Written by Nicola Crosetto in Nov 2025

%[]{}

cd('/Users/nicolacrosetto/Library/CloudStorage/OneDrive-KarolinskaInstitutet/KI/Manuscripts/Li Wang et al_FOLD-FISH/FIRST SUBMISSION/COSMIC-OncoKB probes');

% Load the probe database
filename = 'merged_allprobes_with_cosmicID_final_seq_RS_271025.txt';
data = readtable(filename, 'Delimiter', '\t', 'ReadVariableNames', true);


%% Extract probe features
% List chromosome names in data
chr_list = unique(data.chromosome);

% Extract # of probes per chromosome
for i=1:numel(chr_list)
    probe_per_chr_num(i,1) = numel(unique(data.GeneName(strcmp(data.chromosome, chr_list{i}))));
end

% List all probes
probe_name = unique(data.GeneName);

% For each probe, find # of oligos and probe size
for i=1:numel(probe_name)

    % Index oligos in probe
    oligo_idx = find(strcmp(data.GeneName, probe_name{i}));

    probe_chr{i,1} = data.chromosome{oligo_idx(1)};
    probe_oligo_num{i,1} = numel(oligo_idx);

    % Calculate probe size
    probe_start{i,1} = data.start(oligo_idx(1));
    probe_end{i,1} = data.xEnd(oligo_idx(end));
    probe_size{i,1} = probe_end{i,1} - probe_start{i,1} + 1;

    % Calculate mean, median and std of distance between consecutive oligos
    % in the i-th probe
    inter_oligo_dist = diff(data.start(oligo_idx));
    probe_inter_oligo_mean_dist{i,1} = mean(inter_oligo_dist);
    probe_inter_oligo_median_dist{i,1} = median(inter_oligo_dist);
    probe_inter_oligo_std_dist{i,1} = std(inter_oligo_dist);

    % Calculate mean GC-content of oligos in the i-th probe
    probe_mean_gc{i,1} = mean(data.gc_content(oligo_idx));


end

% Create a table with number of probes per chromosome
probe_per_chr = table(chr_list, probe_per_chr_num, ...
                          'VariableNames', {'Chromosome', 'NumProbes'});

% Merge probe data in a table
probe_data = table(probe_name, probe_chr, probe_oligo_num, probe_start, ...
                   probe_end, probe_size, probe_inter_oligo_mean_dist, ...
                   probe_inter_oligo_median_dist, probe_inter_oligo_std_dist, probe_mean_gc, ...
                   'VariableNames', {'GeneName', 'ProbeChr', 'OligoNum', ...
                                     'Start', 'End', 'Size', ...
                                     'MeanInterOligoDist', 'MedianInterOligoDist', ...
                                     'StdInterOligoDist', 'gc_content'});

% Save probe_data and probe_per_chr_table to a single .mat file
save('probe_data.mat', 'probe_data', 'probe_per_chr');

clear all
clc

