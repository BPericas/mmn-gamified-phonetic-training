%% S04_cluster_permutation_MMN.m
%
%  Fourth step in the MMN analysis pipeline.
%
%  Tests whether the MMN difference wave (DEV - STD) is significantly
%  different from zero across participants using a cluster-based permutation
%  test implemented in FieldTrip (Maris & Oostenveld, 2007). It also compares
%  pre- and post-test MMNs to see if there are significant differences in
%  the MMNs.
%
%  The test compares the real MMN data against a null condition of
%  zero-filled arrays, using a one-tailed dependent-samples t-test
%  (negative tail only, since the MMN is a negative deflection).
%  Cluster correction controls the family-wise error rate across the
%  channel x time space.
%
%  Input:   <subjectID>_ERPdiff.mat files in dataDir
%           (output of S03_extract_ERPs_and_MMN.m)
%  Output:  stat.mat — FieldTrip statistics structure (MMN vs zero, one-tailed)
%           stats_prepostdiff.mat — FieldTrip statistics structure (Pre vs Post,
%                                  two-tailed)
%           Figures (publication-ready format): Sig_clusters.fig (grand-average MMN
%                                               with significant cluster shaded),
%                                               topoplot.fig (topographic t-values),
%                                               PreVsPostDiff.fig (post-minus-pre diff
%                                               wave with significant clusters)
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
dataDir = 'U:\SHaPS\Paul_Iverson\Bego_PhD\gamified_training\ERP_differences';   
% Channel to display in the ERP waveform plot (Plot 1)
% Run  disp(stat.label)  after loading your data to see available labels
plotChannel = 'FCz';

% Epoch time window (seconds)
latencyWindow = [-0.2  0.8];

% Analysis time window (seconds) for the permutation test
tWin = [0  0.8];

% Number of permutations (1000 is fast; 10000 is recommended for publication)
nPerms = 10000;

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
cfg.latency          = tWin;
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
%
%  grandavg spans the full epoch, but stat.time is restricted to cfg.latency
%  (tWin) — so the ERP waveform and the significance shading need their own
%  time vectors rather than sharing one.

erp_time_ms  = grandavg.time * 1000;
stat_time_ms = stat.time * 1000;
sig_times    = any(stat.mask, 1);

ch_idx = find(strcmp(stat.label, plotChannel));
if isempty(ch_idx)
    warning('Channel "%s" not found — using first channel instead. Check plotChannel in USER SETTINGS.', plotChannel);
    ch_idx = 1;
end

% Publication-ready formatting settings (Arial font, standardised sizes)
fontName = 'Arial';
fontSizeLabel = 10;   % axis labels
fontSizeTick = 9;     % axis tick labels

figure('Color', 'white');
hold on;

% Shade significant time windows in grey
ylims   = [min(grandavg.avg(ch_idx,:)) * 1.4, max(grandavg.avg(ch_idx,:)) * 1.4];
sig_on  = find(diff([0, sig_times]) ==  1);
sig_off = find(diff([sig_times,  0]) == -1);
for k = 1:length(sig_on)
    patch([stat_time_ms(sig_on(k))  stat_time_ms(sig_off(k)) ...
           stat_time_ms(sig_off(k)) stat_time_ms(sig_on(k))], ...
          [ylims(1) ylims(1) ylims(2) ylims(2)], ...
          [0.85 0.85 0.85], 'EdgeColor', 'none', 'FaceAlpha', 0.6);
end

plot(erp_time_ms, grandavg.avg(ch_idx,:), 'k', 'LineWidth', 2);
yline(0, '--k', 'LineWidth', 0.8);
xline(0, ':k',  'LineWidth', 0.8);
xlabel('Time (ms)', 'FontName', fontName, 'FontSize', fontSizeLabel, 'FontWeight', 'bold');
ylabel('Amplitude (\muV)', 'FontName', fontName, 'FontSize', fontSizeLabel, 'FontWeight', 'bold');
set(gca, 'FontName', fontName, 'FontSize', fontSizeTick);
ylim(ylims);
xlim([erp_time_ms(1) erp_time_ms(end)]);
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

    figure('Color', 'white');
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
    set(findall(gcf, '-property', 'FontName'), 'FontName', 'Arial');
    set(findall(gcf, '-property', 'FontSize'), 'FontSize', 9);
    savefig(gcf, fullfile(dataDir, 'topoplot.fig'));

else
    fprintf('No significant negative cluster found — topoplot skipped.\n');
end

%% -----------------------------------------------------------------------
%% PLOT 3: Overview of all clusters (ft_clusterplot)

% cfg_cp           = [];
% cfg_cp.alpha     = clAlpha;
% cfg_cp.parameter = 'stat';
% cfg_cp.zlim      = topoZlim;
% ft_clusterplot(cfg_cp, stat);
% savefig(gcf, fullfile(dataDir, 'all_clusts.fig'));

fprintf('\nDone. Figures saved to:\n  %s\n', dataDir);

%% ============================================================
%% ANALYSIS 2: Is post-test MMN amplitude greater than pre-test?
%% ============================================================

% --- Separate pre and post files by naming convention --------
%   Pre-test:  PT*_ERPdiff.mat            (no 'post_test' in name)
%   Post-test: PT*_post_test_ERPdiff.mat
all_files  = dir(fullfile(dataDir, 'PT*_ERPdiff.mat'));
post_files = all_files(contains({all_files.name}, '_post_test'));

% Pre files = all ERPdiff files that do NOT match the post pattern
all_names  = {all_files.name};
post_names = {post_files.name};
pre_files  = all_files(~ismember(all_names, post_names));

% --- Match participants by ID (the PT* prefix before first underscore) ---
pre_ids  = cellfun(@(n) regexp(n, '^(PT[^_]+)', 'match', 'once'), ...
                   {pre_files.name},  'UniformOutput', false);
post_ids = cellfun(@(n) regexp(n, '^(PT[^_]+)', 'match', 'once'), ...
                   {post_files.name}, 'UniformOutput', false);

% --- Keep only participants present in BOTH sessions ----------
[matched_ids, pre_idx, post_idx] = intersect(pre_ids, post_ids, 'stable');
Nsubj_paired = numel(matched_ids);
fprintf('%d participants have both pre- and post-test data.\n', Nsubj_paired);
if Nsubj_paired == 0
    error('No participants with both pre- and post-test ERPdiff files were found.');
end

% --- Load matched pre- and post-test data --------------------
allpre  = cell(1, Nsubj_paired);
allpost = cell(1, Nsubj_paired);

for i = 1:Nsubj_paired
    preFile  = pre_files(pre_idx(i));
    d = load(fullfile(preFile.folder, preFile.name));
    allpre{i}        = [];
    allpre{i}.avg    = d.ERPdiff.erp;
    allpre{i}.time   = d.ERPdiff.times / 1000;
    allpre{i}.label  = {d.ERPdiff.chanlocs.labels}';
    allpre{i}.dimord = 'chan_time';

    postFile = post_files(post_idx(i));
    d = load(fullfile(postFile.folder, postFile.name));
    allpost{i}        = [];
    allpost{i}.avg    = d.ERPdiff.erp;
    allpost{i}.time   = d.ERPdiff.times / 1000;
    allpost{i}.label  = {d.ERPdiff.chanlocs.labels}';
    allpost{i}.dimord = 'chan_time';
end

% --- Configure the test --------------------------------------
% depsamplesT computes t = condition1 - condition2 = post - pre.
% A larger MMN in post means post is more negative, so post - pre < 0
% -> one-tailed test in the negative direction (tail = -1).
cfg2                  = [];
cfg2.channel          = {'all'};
cfg2.latency          = tWin;
cfg2.method           = 'montecarlo';
cfg2.statistic        = 'depsamplesT';
cfg2.correctm         = 'cluster';
cfg2.clusteralpha     = clAlpha;
cfg2.clusterstatistic = 'maxsum';
cfg2.minnbchan        = 2;
cfg2.neighbours       = neighbours;
cfg2.tail             = 0;    % two-tailed: captures both negative (post > pre) and positive (pre > post) clusters
cfg2.clustertail      = 0;
cfg2.alpha            = clAlpha;
cfg2.numrandomization = nPerms;

design2      = zeros(2, Nsubj_paired*2);
design2(1,:) = [1:Nsubj_paired, 1:Nsubj_paired];              % subject IDs (uvar)
design2(2,:) = [ones(1,Nsubj_paired), ones(1,Nsubj_paired)*2]; % cond1=post, cond2=pre (ivar)
cfg2.design  = design2;
cfg2.uvar    = 1;
cfg2.ivar    = 2;

% Run: post (condition 1) vs pre (condition 2)
stat2_prepostdiff      = ft_timelockstatistics(cfg2, allpost{:}, allpre{:});
stat2_prepostdiff.elec = elec;
save(fullfile(dataDir,'stats_prepostdiff.mat'), 'stat2_prepostdiff');
fprintf('Statistics saved to stat2_prepostdiff.mat\n');

%% -----------------------------------------------------------
%% PLOT: Post-minus-pre difference wave with significant clusters
%% -----------------------------------------------------------

cfg_ga2 = [];
cfg_ga2.keepindividual = 'no';

grandavg_post = ft_timelockgrandaverage(cfg_ga2, allpost{:});
grandavg_pre = ft_timelockgrandaverage(cfg_ga2, allpre{:});

diff_avg = grandavg_post.avg - grandavg_pre.avg; % post - pre

ch_name = plotChannel;

% Time vectors
erp_time_ms = grandavg_post.time * 1000;
stat_time_ms = stat2_prepostdiff.time * 1000;

% Find requested channel separately in ERP and stat structures
ch_idx_erp = find(strcmp(grandavg_post.label, ch_name), 1);
ch_idx_stat = find(strcmp(stat2_prepostdiff.label, ch_name), 1);

if isempty(ch_idx_erp)
error('Channel %s not found in grand-average ERP data.', ch_name);
end

if isempty(ch_idx_stat)
error('Channel %s not found in statistical result.', ch_name);
end

% Significant time points across the statistical analysis
sig_times2 = any(stat2_prepostdiff.mask, 1);

% Publication-ready formatting: Arial font with standardised sizes
fontName = 'Arial';
fontSizeLabel = 10;   % axis labels
fontSizeTick = 9;     % axis tick labels

figure('Color', 'white');
hold on;

wave = diff_avg(ch_idx_erp,:);

ylims2 = [min(wave)*1.4, max(wave)*1.4];

% Identify continuous significant intervals
sig_on2 = find(diff([0, sig_times2]) == 1);
sig_off2 = find(diff([sig_times2, 0]) == -1);

% Shade significant windows using the STATISTICAL time vector
for k = 1:length(sig_on2)
    t1 = stat_time_ms(sig_on2(k));
    t2 = stat_time_ms(sig_off2(k));
    patch([t1 t2 t2 t1], ...
        [ylims2(1) ylims2(1) ylims2(2) ylims2(2)], ...
        [0.85 0.85 0.85], ...
        'EdgeColor', 'none', ...
        'FaceAlpha', 0.6);
end

% Plot ERP using the ERP time vector
plot(erp_time_ms, wave, 'b', 'LineWidth', 2);
yline(0, '--k', 'LineWidth', 0.8);
xline(0, ':k', 'LineWidth', 0.8);
xlabel('Time (ms)', 'FontName', fontName, 'FontSize', fontSizeLabel, 'FontWeight', 'bold');
ylabel('\DeltaAmplitude (\muV) [Post - Pre]', 'FontName', fontName, 'FontSize', fontSizeLabel, 'FontWeight', 'bold');
set(gca, 'FontName', fontName, 'FontSize', fontSizeTick, 'YDir', 'reverse');
ylim(ylims2);
xlim([erp_time_ms(1) erp_time_ms(end)]);
box off;
savefig(gcf, fullfile(dataDir, 'PreVsPostDiff.fig'));