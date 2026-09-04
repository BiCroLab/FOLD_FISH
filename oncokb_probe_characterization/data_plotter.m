%A script to plot the features of FOLD-FISH probes targeting OncoKB and COSMIC genes
%designed by Raquel Silva for the FOLD-FISH manuscript
%
%Written by Nicola Crosetto in Nov 2025

%[]{}

cd('/Users/nicolacrosetto/Library/CloudStorage/OneDrive-KarolinskaInstitutet/KI/Manuscripts/Li Wang et al_FOLD-FISH/FIRST SUBMISSION/COSMIC-OncoKB probes');

% Load the parsed probe data
load('probe_data.mat');

% Plot number of probes per chromosome
figure(1);
bar(probe_per_chr.NumProbes);
xlabel('Chromosome');
ylabel('Number of Probes');


% Plot distribution of number of oligos per probe as a boxplot
figure(2);
boxplot(cell2mat(probe_data.OligoNum), 'Labels', {'Oligos per Probe'}, 'Positions', 1, 'BoxStyle', 'outline', 'Color', 'k');
hold on; % Hold the current plot to overlay data points
% Show data points on top of the boxplot with horizontal jitter
scatter(ones(size(cell2mat(probe_data.OligoNum))) + 0.03 * randn(size(cell2mat(probe_data.OligoNum))), ...
        cell2mat(probe_data.OligoNum), 20, [0.5 0.5 0.5], 'filled', 'MarkerFaceAlpha', 0.5); % Smaller filled gray points
hold off; % Release the hold on the current plot
xlabel('Oligos per Probe');
ylabel('Frequency');
% Calculate and display mean and median
meanOligos = mean(cell2mat(probe_data.OligoNum));
medianOligos = median(cell2mat(probe_data.OligoNum));
annotation('textbox', [0.7, 0.7, 0.2, 0.1], 'String', {sprintf('Mean: %.2f', meanOligos), sprintf('Median: %.2f', medianOligos)}, ...
           'Color', 'k', 'FontSize', 10, 'FitBoxToText', 'on', 'EdgeColor', 'none');


% Plot distribution of mean inter-oligo distance per probe as a histogram
figure(3);
histogram(cell2mat(probe_data.MeanInterOligoDist), 'Normalization', 'probability', 'FaceColor', [0.5 0.5 0.5]);
xlabel('Mean Inter-Oligo Distance');
ylabel('Probability');
% Calculate and display mean and median
meanInterOligoDist = mean(cell2mat(probe_data.MeanInterOligoDist));
medianInterOligoDist = median(cell2mat(probe_data.MeanInterOligoDist));
annotation('textbox', [0.7, 0.7, 0.2, 0.1], 'String', {sprintf('Mean: %.2f', meanInterOligoDist), sprintf('Median: %.2f', medianInterOligoDist)}, ...
           'Color', 'k', 'FontSize', 10, 'FitBoxToText', 'on', 'EdgeColor', 'none');


% Plot distribution of standard deviation of inter-oligo distance per probe as a histogram
figure(4);
histogram(cell2mat(probe_data.StdInterOligoDist), 'Normalization', 'probability', 'FaceColor', [0.5 0.5 0.5]);
xlabel('Standard Deviation of Inter-Oligo Distance');
ylabel('Probability');
% Calculate and display mean and median
meanMeanStd = mean(cell2mat(probe_data.StdInterOligoDist));
medianMeanStd = median(cell2mat(probe_data.StdInterOligoDist));
annotation('textbox', [0.7, 0.7, 0.2, 0.1], 'String', {sprintf('Mean: %.2f', meanMeanStd), sprintf('Median: %.2f', medianMeanStd)}, ...
           'Color', 'k', 'FontSize', 10, 'FitBoxToText', 'on', 'EdgeColor', 'none');


% Plot distribution of probe sizes as a boxplot
figure(5);
boxplot(cell2mat(probe_data.Size), 'Labels', {'Probe Size'}, 'Positions', 1, 'BoxStyle', 'outline', 'Color', 'k');
hold on; % Hold the current plot to overlay data points
% Show data points on top of the boxplot with horizontal jitter
scatter(ones(size(cell2mat(probe_data.Size))) + 0.03 * randn(size(cell2mat(probe_data.Size))), ...
        cell2mat(probe_data.Size), 20, [0.5 0.5 0.5], 'filled', 'MarkerFaceAlpha', 0.5); % Smaller filled gray points
hold off; % Release the hold on the current plot
xlabel('Probe Size');
ylabel('Frequency');
% Calculate and display mean and median
meanSize = mean(cell2mat(probe_data.Size));
medianSize = median(cell2mat(probe_data.Size));
annotation('textbox', [0.7, 0.7, 0.2, 0.1], 'String', {sprintf('Mean: %.2f', meanSize), sprintf('Median: %.2f', medianSize)}, ...
           'Color', 'k', 'FontSize', 10, 'FitBoxToText', 'on', 'EdgeColor', 'none');


% Plot distribution of GC content as a histogram
figure(6);
histogram(cell2mat(probe_data.gc_content), 'Normalization', 'probability', 'FaceColor', [0.5 0.5 0.5]);
xlabel('Mean GC-content per probe');
ylabel('Probability');
% Calculate and display mean and median
meanGCContent = mean(cell2mat(probe_data.gc_content));
medianGCContent = median(cell2mat(probe_data.gc_content));
annotation('textbox', [0.7, 0.7, 0.2, 0.1], 'String', {sprintf('Mean: %.2f', meanGCContent), sprintf('Median: %.2f', medianGCContent)}, ...
           'Color', 'k', 'FontSize', 10, 'FitBoxToText', 'on', 'EdgeColor', 'none');

% Save figures as PDF files with variable names in the filenames
variableNames = {'oncokb_probe_per_chr', 'oncokb_probe_oligo_num', 'oncokb_probe_mean_interoligo_dist', ...
                'oncokb_mean_inter_oligo_std', 'oncokb_probe_size', 'oncokb_probe_gc_content'};
for figNum = 1:6
    saveas(figure(figNum), sprintf('%s.pdf', variableNames{figNum}));
end

close all
