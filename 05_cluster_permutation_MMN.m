%% 05_cluster_permutation_MMN.m
%
%  Fifth step in the MMN analysis pipeline.
%
%  Tests whether the MMN difference wave (DEV - STD) is significantly
%  different from zero across participants using a cluster-based permutation
%  test implemented in FieldTrip (Maris & Oostenveld, 2007).
%
%  The test compares the real MMN data against a null condition of
%  zero-filled arrays, using a one-tailed dependent-samples t-test
%  (negative tail only, since the MMN is a negative deflection).
%  Cluster correction controls the family-wise error rate across the
%  channel x time space.
%
%  Input:   <subjectID>_ERPdiff.mat files in dataDir
%           (output of 03_extract_ERPs_and_MMN.m)
%  Output:  stat.mat — FieldTrip statistics structure
%           Sig_clusters.fig, topoplot.fig, all_clusts.fig — saved figures
%
%  Dependencies
%  ------------
%  - FieldTrip (https://www.fieldtriptoolbox.org)
%
%  Reference
%  ---------
%  Maris, E., & Oostenveld, R. (2007). Nonparametric statistical testing of
%  EEG- and MEG-data. Journal of Neuroscience Methods, 164(1), 177-190.

clear; ft_defaults;

%% -----------------------------------------------------------------------
%% USER SETTINGS — edit these before running

% Folder containing the _ERPdiff.mat files
dataDir = pwd;   % defaults to current folder; replace with an explicit path if needed
% Example: dataDir = 'C:\MyProject\ERP_differences';

% Channel to display in the ERP waveform plot (Plot 1)
% Run  disp(stat.label)  after loading your data to see available labels
plotChannel = 'FCz';

% Analysis time window (seconds)
latencyWindow = [-0.2  0.8];

% Number of permutations (1000 is fast; 10000 is recommended for publication)
nPerms = 1000;

% Cluster-forming and cluster-level significance threshold
clAlpha = 0.05;

% Minimum number of neighbouring channels required to form a cluster
minNbChan = 2;

% Colour scale limit for the cluster overview plot (Plot 3)
% Increase if your t-values exceed this range
topoZlim = [-4  4];

%% -----------------------------------------------------------------------
%% Load subject data and build FieldTrip timelock structures

files  = dir(fullfile(dataDir, '*_ERPdiff.mat'));
nFiles = length(files);

if nFiles == 0
    error('No *_ERPdiff.mat files found in:\n  %s\nCheck dataDir.', dataDir);
end

fprintf('Found %d participant file(s).\n', nFiles);

allsubj = cell(1, nFiles);

for i = 1:nFiles
    d = load(fullfile(files(i).folder, files(i).name));

    ft_data        = [];
    ft_data.avg    = d.ERPdiff.erp;
    ft_data.time   = d.ERPdiff.times / 1000;   % ms to seconds
    ft_data.label  = {d.ERPdiff.chanlocs.labels}';
    ft_data.dimord = 'chan_time';

    allsubj{i} = ft_data;
end

%% -----------------------------------------------------------------------
%% Build electrode neighbour structure from actual electrode positions
%
%  Using the real electrode coordinates from the data rather than a
%  template ensures neighbours are correct for your specific cap layout.

d    = load(fullfile(files(1).folder, files(1).name));
elec = [];
elec.label   = {d.ERPdiff.chanlocs.labels}';
elec.elecpos = [-[d.ERPdiff.chanlocs.Y]', [d.ERPdiff.chanlocs.X]', [d.ERPdiff.chanlocs.Z]'];
elec.chanpos = elec.elecpos;
elec.unit    = 'mm';

cfg_neigh        = [];
cfg_neigh.method = 'triangulation';
cfg_neigh.elec   = elec;
neighbours       = ft_prepare_neighbours(cfg_neigh);

%% -----------------------------------------------------------------------
%% Cluster-based permutation test (MMN vs zero)
%
%  The null condition is a set of zero-filled arrays with the same
%  dimensions as the real data. A significant negative cluster indicates
%  that the MMN is reliably different from zero across participants.

Nsubj = nFiles;

cfg                  = [];
cfg.channel          = {'all'};
cfg.latency          = latencyWindow;
cfg.method           = 'montecarlo';
cfg.statistic        = 'depsamplesT';
cfg.correctm         = 'cluster';
cfg.clusteralpha     = clAlpha;
cfg.clusterstatistic = 'maxsum';
cfg.minnbchan        = minNbChan;
cfg.neighbours       = neighbours;
cfg.tail             = -1;       % one-tailed: negative direction (MMN is negative)
cfg.clustertail      = -1;
cfg.alpha            = clAlpha;
cfg.numrandomization = nPerms;

% Design matrix: row 1 = subject ID, row 2 = condition (1=data, 2=zeros)
design      = zeros(2, Nsubj * 2);
design(1,:) = [1:Nsubj, 1:Nsubj];
design(2,:) = [ones(1,Nsubj), ones(1,Nsubj)*2];
cfg.design  = design;
cfg.uvar    = 1;   % unit variable (subjects)
cfg.ivar    = 2;   % independent variable (condition)

% Create zero-filled null condition
allzeros = allsubj;
for i = 1:Nsubj
    allzeros{i}.avg = zeros(size(allsubj{i}.avg));
end

fprintf('Running cluster permutation test (%d randomisations)...\n', nPerms);
stat = ft_timelockstatistics(cfg, allsubj{:}, allzeros{:});

% Save statistics
save(fullfile(dataDir, 'stat.mat'), 'stat');
fprintf('Statistics saved to stat.mat\n');

%% -----------------------------------------------------------------------
%% Add electrode positions to stat (required for topographic plots)
stat.elec = elec;

%% -----------------------------------------------------------------------
%% Compute grand-average difference wave across all participants

cfg_ga              = [];
cfg_ga.keepindividual = 'no';
grandavg = ft_timelockgrandaverage(cfg_ga, allsubj{:});

%% -----------------------------------------------------------------------
%% PLOT 1: Grand-average ERP waveform with significant cluster(s) shaded

time_ms   = stat.time * 1000;
sig_times = any(stat.mask, 1);

ch_idx = find(strcmp(stat.label, plotChannel));
if isempty(ch_idx)
    warning('Channel "%s" not found — using first channel instead. Check plotChannel in USER SETTINGS.', plotChannel);
    ch_idx = 1;
end

figure('Name', 'MMN Difference Wave');
hold on;

% Shade significant time windows in grey
ylims   = [min(grandavg.avg(ch_idx,:)) * 1.4, max(grandavg.avg(ch_idx,:)) * 1.4];
sig_on  = find(diff([0, sig_times]) ==  1);
sig_off = find(diff([sig_times,  0]) == -1);
for k = 1:length(sig_on)
    patch([time_ms(sig_on(k))  time_ms(sig_off(k)) ...
           time_ms(sig_off(k)) time_ms(sig_on(k))], ...
          [ylims(1) ylims(1) ylims(2) ylims(2)], ...
          [0.85 0.85 0.85], 'EdgeColor', 'none', 'FaceAlpha', 0.6);
end

plot(time_ms, grandavg.avg(ch_idx,:), 'k', 'LineWidth', 2);
yline(0, '--k', 'LineWidth', 0.8);
xline(0, ':k',  'LineWidth', 0.8);
xlabel('Time (ms)');
ylabel('Amplitude (\muV)');
title(['MMN at ' plotChannel ' — shaded = significant cluster(s)']);
ylim(ylims);
xlim([time_ms(1) time_ms(end)]);
set(gca, 'YDir', 'reverse');   % ERP convention: negative up
box off;
savefig(gcf, fullfile(dataDir, 'Sig_clusters.fig'));

%% -----------------------------------------------------------------------
%% PLOT 2: Topographic map of t-values for the main significant cluster
%
%  Only runs if at least one significant negative cluster was found.

if isfield(stat, 'negclusters') && ~isempty(stat.negclusters) && ...
   stat.negclusters(1).prob < clAlpha

    c1_time_idx = any(stat.negclusterslabelmat == 1, 1);
    c1_window   = [stat.time(find(c1_time_idx,1,'first')), ...
                   stat.time(find(c1_time_idx,1,'last'))];

    figure('Name', 'Topomap — Cluster 1');
    cfg_topo                  = [];
    cfg_topo.parameter        = 'stat';
    cfg_topo.xlim             = c1_window;
    cfg_topo.zlim             = 'maxabs';
    cfg_topo.colormap         = 'RdBu';
    cfg_topo.highlight        = 'on';
    cfg_topo.highlightchannel = stat.label(any(stat.negclusterslabelmat == 1, 2));
    cfg_topo.highlightsymbol  = '.';
    cfg_topo.highlightsize    = 14;
    cfg_topo.colorbar         = 'yes';
    cfg_topo.comment          = ['p = ' num2str(stat.negclusters(1).prob, '%.4f')];
    cfg_topo.commentpos       = 'rightbottom';
    ft_topoplotER(cfg_topo, stat);
    title(['T-values ' num2str(c1_window(1)*1000,'%.0f') '-' ...
           num2str(c1_window(2)*1000,'%.0f') ' ms | dots = cluster channels']);
    savefig(gcf, fullfile(dataDir, 'topoplot.fig'));

else
    fprintf('No significant negative cluster found — topoplot skipped.\n');
end

%% -----------------------------------------------------------------------
%% PLOT 3: Overview of all clusters (ft_clusterplot)

cfg_cp           = [];
cfg_cp.alpha     = clAlpha;
cfg_cp.parameter = 'stat';
cfg_cp.zlim      = topoZlim;
ft_clusterplot(cfg_cp, stat);
savefig(gcf, fullfile(dataDir, 'all_clusts.fig'));

fprintf('\nDone. Figures saved to:\n  %s\n', dataDir);
