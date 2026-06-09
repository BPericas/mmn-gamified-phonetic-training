%% 06_cluster_permutation_pre_vs_post.m
%
%  Sixth and final step in the MMN analysis pipeline.
%
%  Tests whether the MMN changed significantly from pre-test to post-test,
%  separately for each phonetic contrast, using a cluster-based permutation
%  test implemented in FieldTrip (Maris & Oostenveld, 2007).
%
%  A two-tailed dependent-samples t-test is run over a 2-D channel x time
%  space. Only subjects with data for both sessions are included.
%
%  Sign convention
%  ---------------
%  depsamplesT computes Post - Pre. Because the MMN is a negative
%  deflection, a training effect (larger MMN after training) produces
%  Post < Pre, yielding a NEGATIVE t-value and a NEGATIVE cluster.
%    Negative cluster = MMN larger at post-test  = improvement
%    Positive cluster = MMN smaller at post-test = no improvement / regression
%
%  Input:   <subjectID>_pre_MMN_<contrast>.mat  \
%           <subjectID>_post_MMN_<contrast>.mat  > from 04_MMN_by_contrast_pre_post.m
%  Output:  ClusterPerm_PreVsPost_FT.mat — FieldTrip stat structs, one per contrast
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

% Folder containing the pre/post MMN .mat files (output of script 04)
outDir = 'C:\MyProject\ERP_differences_by_contrast';   % <-- change to your folder

% Contrast identifiers — must match the suffixes in your filenames
% e.g. files named sub01_pre_MMN_AB.mat use contrast label 'AB'
contrasts = {'AB', 'CD'};

% Human-readable contrast labels used in plot titles
% Adapt to your phonetic contrasts
contrastLabels = containers.Map({'AB','CD'}, {'/i:/ - /I/', '/ae/ - /e/'});

% Frontocentral channels to include in the analysis (FieldTrip mixed-case convention)
% Adjust to match your electrode layout
chanLabels = {'F3','F4','FC3','FCz','FC4','C3','Cz','C4'};

% Analysis time window (seconds) — adjust after inspecting grand averages
tWin = [0  0.5];

% Number of permutations (10000 recommended for publication)
nPerms = 10000;

% Cluster-forming and cluster-level significance threshold
clAlpha = 0.05;

% Colour scheme for grand-average plots (one colour per contrast)
% Add entries if you have more than two contrasts
contrastColours = containers.Map({'AB','CD'}, {[0 0 0.8], [0.8 0 0]});

%% -----------------------------------------------------------------------
%% Build channel neighbour structure from actual electrode positions
%
%  Using real electrode coordinates rather than a template ensures
%  neighbours are correct for your specific cap layout.

% Find a representative file to extract electrode positions
repList = [];
for ci = 1:numel(contrasts)
    repList = dir(fullfile(outDir, sprintf('*_pre_MMN_%s.mat', contrasts{ci})));
    if ~isempty(repList), break; end
end
if isempty(repList)
    error('No pre-test MMN files found in:\n  %s\nRun script 04 first.', outDir);
end

tmp   = load(fullfile(outDir, repList(1).name), 'd');
d_rep = tmp.d;

xyz = [[d_rep.chanlocs.X]', [d_rep.chanlocs.Y]', [d_rep.chanlocs.Z]'];
if any(isnan(xyz(:))) || all(xyz(:) == 0)
    warning('Electrode XYZ positions missing — falling back to EEG1010 template.');
    cfg_nb          = [];
    cfg_nb.method   = 'template';
    cfg_nb.template = 'EEG1010_neighb.mat';
else
    elec.label   = normaliseLabels({d_rep.chanlocs.labels}');
    elec.chanpos = xyz;
    elec.elecpos = xyz;
    elec.unit    = 'mm';
    cfg_nb        = [];
    cfg_nb.method = 'distance';
    cfg_nb.elec   = elec;
end
neighbours = ft_prepare_neighbours(cfg_nb);
fprintf('Neighbours computed for %d channels.\n', numel(d_rep.chanlocs));

%% -----------------------------------------------------------------------
%% Main loop: one cluster permutation test per contrast

results = struct();

for ci = 1:numel(contrasts)

    cont      = contrasts{ci};
    contLabel = contrastLabels(cont);
    fprintf('\n======= CONTRAST: %s (%s) =======\n', cont, contLabel);

    %% Find subjects with BOTH pre and post files
    preFiles  = dir(fullfile(outDir, sprintf('*_pre_MMN_%s.mat',  cont)));
    postFiles = dir(fullfile(outDir, sprintf('*_post_MMN_%s.mat', cont)));

    % Extract subject IDs from filenames
    preIDs  = cellfun(@(n) regexp(n, '^(.+?)_pre_MMN',  'tokens', 'once'), {preFiles.name},  'UniformOutput', false);
    postIDs = cellfun(@(n) regexp(n, '^(.+?)_post_MMN', 'tokens', 'once'), {postFiles.name}, 'UniformOutput', false);
    preIDs  = [preIDs{:}];
    postIDs = [postIDs{:}];

    pairedIDs = intersect(preIDs, postIDs, 'stable');
    nPaired   = numel(pairedIDs);
    fprintf('  %d subjects with both pre and post data.\n', nPaired);

    if nPaired < 3
        warning('Too few paired subjects for contrast %s — skipping.', cont);
        continue;
    end

    %% Load paired FieldTrip timelocks
    ftPre  = {};
    ftPost = {};

    for s = 1:nPaired
        fPre  = fullfile(outDir, sprintf('%s_pre_MMN_%s.mat',  pairedIDs{s}, cont));
        fPost = fullfile(outDir, sprintf('%s_post_MMN_%s.mat', pairedIDs{s}, cont));

        tmp = load(fPre,  'd');  ftPre{end+1}  = toFieldTrip(tmp.d);  %#ok<SAGROW>
        tmp = load(fPost, 'd');  ftPost{end+1} = toFieldTrip(tmp.d);  %#ok<SAGROW>
    end

    %% Run cluster-based permutation test (two-tailed, pre vs post)
    cfg                  = [];
    cfg.method           = 'montecarlo';
    cfg.statistic        = 'depsamplesT';   % paired t-test
    cfg.correctm         = 'cluster';
    cfg.clusteralpha     = clAlpha;
    cfg.clusterstatistic = 'maxsum';
    cfg.tail             = 0;               % two-tailed
    cfg.clustertail      = 0;
    cfg.correcttail      = 'prob';          % correct p-values for two-tailed test
    cfg.alpha            = clAlpha;
    cfg.numrandomization = nPerms;
    cfg.neighbours       = neighbours;
    cfg.channel          = chanLabels;
    cfg.latency          = tWin;

    % Design matrix:
    %   Row 1 (ivar) — session:  1 = pre,  2 = post
    %   Row 2 (uvar) — subject:  1..nPaired, repeated for each session
    % Data order must match: all pre subjects first, then all post
    cfg.design = [ones(1,nPaired), 2*ones(1,nPaired); ...
                  1:nPaired,       1:nPaired         ];
    cfg.ivar   = 1;
    cfg.uvar   = 2;

    fprintf('  Running permutation test (%d randomisations)...\n', nPerms);
    stat = ft_timelockstatistics(cfg, ftPre{:}, ftPost{:});

    %% Print cluster summary to console
    printClusters(stat, contLabel, nPaired, clAlpha);

    %% Generate plots
    plotTStatistic(stat, contLabel, nPaired, clAlpha);
    plotGrandAverage(ftPre, ftPost, stat, cont, contLabel, nPaired, chanLabels, clAlpha, contrastColours);

    %% Store results
    results.(cont)               = stat;
    results.([cont '_subjects']) = pairedIDs;

end  % contrast loop

save(fullfile(outDir, 'ClusterPerm_PreVsPost_FT.mat'), 'results');
fprintf('\nResults saved to:\n  %s\n', fullfile(outDir, 'ClusterPerm_PreVsPost_FT.mat'));

%% -----------------------------------------------------------------------
%% LOCAL FUNCTIONS
%% -----------------------------------------------------------------------

function tl = toFieldTrip(d)
% Convert an EEGLAB-style MMN struct to a FieldTrip timelock structure.
    tl.avg    = d.erp;
    tl.time   = d.times / 1000;   % ms to seconds
    tl.label  = normaliseLabels({d.chanlocs.labels}');
    tl.dimord = 'chan_time';
end

% -------------------------------------------------------------------------
function labels = normaliseLabels(labels)
% Convert all-caps z-suffix labels to mixed-case (e.g. FCZ -> FCz, CZ -> Cz).
% Required because EEGLAB/NeuroScan exports all-caps while FieldTrip uses
% mixed-case. Remove or adapt if your system already uses mixed-case labels.
    for k = 1:numel(labels)
        lbl = labels{k};
        if numel(lbl) >= 2 && lbl(end) == 'Z'
            lbl(end)  = 'z';
            labels{k} = lbl;
        end
    end
end

% -------------------------------------------------------------------------
function printClusters(stat, contLabel, nPaired, clAlpha)
% Print a summary of significant clusters to the console.
    fprintf('\n  ---- Pre vs Post: %s  (N = %d paired subjects) ----\n', contLabel, nPaired);
    fprintf('  Negative clusters (Post < Pre amplitude = MMN LARGER at post = improvement):\n');
    if isfield(stat,'negclusters') && ~isempty(stat.negclusters)
        for k = 1:numel(stat.negclusters)
            p   = stat.negclusters(k).prob;
            msk = stat.negclusterslabelmat == k;
            t1  = stat.time(find(any(msk,1),1,'first')) * 1000;
            t2  = stat.time(find(any(msk,1),1,'last'))  * 1000;
            sig = ''; if p < clAlpha, sig = '  ***'; end
            fprintf('    Cluster %d:  p = %.3f%s,  [%.0f - %.0f ms]\n', k, p, sig, t1, t2);
        end
    else
        fprintf('    None above threshold.\n');
    end
    fprintf('  Positive clusters (Post > Pre amplitude = MMN SMALLER at post = no improvement):\n');
    if isfield(stat,'posclusters') && ~isempty(stat.posclusters)
        for k = 1:numel(stat.posclusters)
            p   = stat.posclusters(k).prob;
            msk = stat.posclusterslabelmat == k;
            t1  = stat.time(find(any(msk,1),1,'first')) * 1000;
            t2  = stat.time(find(any(msk,1),1,'last'))  * 1000;
            sig = ''; if p < clAlpha, sig = '  ***'; end
            fprintf('    Cluster %d:  p = %.3f%s,  [%.0f - %.0f ms]\n', k, p, sig, t1, t2);
        end
    else
        fprintf('    None above threshold.\n');
    end
    fprintf('\n');
end

% -------------------------------------------------------------------------
function plotTStatistic(stat, contLabel, nPaired, clAlpha)
% Plot channel-averaged t-statistic over time with significant cluster
% windows shaded. Green = negative cluster (improvement), orange = positive.
    tMean = mean(stat.stat, 1);
    tMs   = stat.time * 1000;

    % Cluster-forming t-threshold (requires Statistics Toolbox; falls back to 2)
    try
        tCrit = tinv(1 - clAlpha/2, nPaired - 1);
    catch
        tCrit = 2;
        warning('tinv not available — using tCrit = 2 as fallback.');
    end
    yLim = max(max(abs(tMean)) * 1.4, tCrit * 1.4);

    figure('Name', sprintf('t-statistic: Pre vs Post — %s  (N = %d)', contLabel, nPaired));
    hold on;

    if isfield(stat,'negclusters')
        for k = 1:numel(stat.negclusters)
            if stat.negclusters(k).prob < clAlpha
                shadeBand(tMs(any(stat.negclusterslabelmat == k, 1)), yLim, [0 0.6 0]);
            end
        end
    end
    if isfield(stat,'posclusters')
        for k = 1:numel(stat.posclusters)
            if stat.posclusters(k).prob < clAlpha
                shadeBand(tMs(any(stat.posclusterslabelmat == k, 1)), yLim, [1 0.5 0]);
            end
        end
    end

    plot(tMs, tMean, 'k-', 'LineWidth', 2);
    yline( tCrit, 'k--', 'LineWidth', 1, 'Label', sprintf('p<%.2f', clAlpha));
    yline(-tCrit, 'k--', 'LineWidth', 1);
    yline(0, 'Color', [0.5 0.5 0.5], 'LineStyle', ':');
    xline(0, 'Color', [0.5 0.5 0.5], 'LineStyle', ':');
    ylim([-yLim yLim]); xlim([tMs(1) tMs(end)]);
    xlabel('Time (ms)'); ylabel('t-statistic (channel average)');
    title(sprintf('Pre vs Post — %s  (N = %d paired subjects)', contLabel, nPaired));
    grid on; hold off;
end

% -------------------------------------------------------------------------
function plotGrandAverage(ftPre, ftPost, stat, cont, contLabel, nPaired, chanLabels, clAlpha, colours)
% Grand-average MMN waveform for pre and post with cluster windows shaded.
    nCh = size(ftPre{1}.avg, 1);
    nT  = size(ftPre{1}.avg, 2);
    stackPre  = zeros(nCh, nT, nPaired);
    stackPost = zeros(nCh, nT, nPaired);
    for s = 1:nPaired
        stackPre(:,:,s)  = ftPre{s}.avg;
        stackPost(:,:,s) = ftPost{s}.avg;
    end
    gaPre  = mean(stackPre,  3);
    gaPost = mean(stackPost, 3);

    times  = ftPre{1}.time * 1000;
    fcIdx  = find(ismember(upper(ftPre{1}.label), upper(chanLabels)));
    preMean  = mean(gaPre(fcIdx,  :), 1);
    postMean = mean(gaPost(fcIdx, :), 1);

    baseCol = colours(cont);
    yLim    = max(max(abs([preMean, postMean]))) * 1.4;
    if yLim == 0, yLim = 1; end

    figure('Name', sprintf('Grand Average MMN — %s (N = %d)', contLabel, nPaired));
    hold on;

    if isfield(stat,'negclusters')
        for k = 1:numel(stat.negclusters)
            if stat.negclusters(k).prob < clAlpha
                shadeBand(stat.time(any(stat.negclusterslabelmat==k,1))*1000, yLim, [0 0.6 0]);
            end
        end
    end
    if isfield(stat,'posclusters')
        for k = 1:numel(stat.posclusters)
            if stat.posclusters(k).prob < clAlpha
                shadeBand(stat.time(any(stat.posclusterslabelmat==k,1))*1000, yLim, [1 0.5 0]);
            end
        end
    end

    hPre  = plot(times, preMean,  '-',  'Color', baseCol, 'LineWidth', 2);
    hPost = plot(times, postMean, '--', 'Color', baseCol, 'LineWidth', 2);
    xline(0,'k:','LineWidth',1); yline(0,'k:','LineWidth',1);
    xlim([times(1) times(end)]); ylim([-yLim yLim]);
    xlabel('Time (ms)'); ylabel('Amplitude (\muV)');
    legend([hPre, hPost], {sprintf('Pre  (n=%d)', nPaired), sprintf('Post (n=%d)', nPaired)}, 'Location','best');
    title(sprintf('Grand Average MMN — %s  (N = %d paired subjects)', contLabel, nPaired));
    grid on; hold off;
end

% -------------------------------------------------------------------------
function shadeBand(tMs, yLim, colour)
% Shade a time band with the given colour.
    if isempty(tMs), return; end
    patch([tMs(1) tMs(end) tMs(end) tMs(1)], [-yLim -yLim yLim yLim], ...
          colour, 'FaceAlpha', 0.20, 'EdgeColor', 'none');
end
