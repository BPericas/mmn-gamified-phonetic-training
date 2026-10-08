%% S05_cluster_permutation_within_phoneme.m
%
%  Fifth step in the MMN analysis pipeline — identity/acoustic confound
%  control analysis.
%
%  For each individual vowel (bit, beat, bet, bat), this script compares
%  the ERP evoked by that vowel when it played the STANDARD role against
%  the ERP evoked by the very same vowel when it played the DEVIANT role
%  (in a different block). Because the physical stimulus is identical in
%  both cases, a genuine difference confirms that the MMN reflects a
%  context/prediction-error effect rather than a low-level acoustic
%  difference between the two vowels in a contrast.
%
%  This mirrors S06_cluster_permutation_within_contrast.m, but the paired
%  comparison is STD-role vs DEV-role for one vowel instead of pre vs post
%  for one contrast.
%
%  Block design
%  ------------
%  Each session (pre/post) consists of 4 blocks, A-D. In each block one
%  vowel is the standard (marker 1) and the other is the deviant (marker 2):
%    Block A: bit  = standard, beat = deviant
%    Block B: beat = standard, bit  = deviant
%    Block C: bet  = standard, bat  = deviant
%    Block D: bat  = standard, bet  = deviant
%  Change stdPhonemeOfBlock / devPhonemeOfBlock in USER SETTINGS if your
%  block-to-role assignment differs.
%
%  Sign convention
%  ---------------
%  depsamplesT computes DEViant-role minus STandard-role (condition 2 minus
%  condition 1).
%    Negative cluster = vowel sounds MORE negative when it is the deviant
%                        = genuine context effect (replicates the MMN
%                          pattern even though the physical stimulus is
%                          unchanged)
%    Positive cluster = vowel sounds LESS negative / more positive when
%                        it is the deviant = unexpected direction
%
%  Only subjects who completed both pre-test and post-test are included, 
%  and the test is run separately for the pre-test session and the post-test session.
%
%  A second analysis then builds the MMN difference wave (DEV-role minus
%  STD-role) for each phoneme and session, and runs a Pre vs Post cluster
%  permutation test on that difference wave (same logic as
%  S06_cluster_permutation_within_contrast.m, but per phoneme instead of
%  per contrast). depsamplesT computes Post - Pre, so a negative cluster
%  means the phoneme's MMN got larger (more negative) after training.
%
%
%  Input:   <subjectID>_GEDAIclean.set / <subjectID>_post..._GEDAIclean.set
%           pairblock_summary.csv / pairblock_summary_post.csv
%  Output:  ERP_differences_by_phoneme/<subjectID>_<session>_role_<phoneme>.mat
%           ClusterPerm_StdVsDev_byPhoneme_FT.mat — FieldTrip stat structs,
%                                                    one per session x phoneme
%           ClusterPerm_PreVsPost_byPhoneme_FT.mat — FieldTrip stat structs,
%                                                     one per phoneme
%           Figures (publication-ready format): per-phoneme grand-average ERPs
%                                               and t-statistic plots with cluster
%                                               overlays; pre vs post MMN variants
%
%  Dependencies
%  ------------
%  - FieldTrip (https://www.fieldtriptoolbox.org)
%  - EEGLAB
%
%  Reference
%  ---------
%  Maris, E., & Oostenveld, R. (2007). Nonparametric statistical testing of
%  EEG- and MEG-data. Journal of Neuroscience Methods, 164(1), 177-190.

clear; eeglab; ft_defaults;

%% -----------------------------------------------------------------------
%% USER SETTINGS — edit these before running

% Root directory of the project (contains the _GEDAIclean.set files)
baseDir = pwd;   % <-- change to your project folder

% File patterns to search for in each session folder
filePatterns = {'*_GEDAIclean.set'};

% CSV files mapping each participant to their block order
% Each file must have columns: Participant_ID, Block_1, Block_2, Block_3, Block_4
preBlockTable  = fullfile(baseDir, 'pairblock_summary.csv');
postBlockTable = fullfile(baseDir, 'pairblock_summary_post.csv');

% Output directory for subject role files and cluster stats
outDir = fullfile(baseDir, 'ERP_differences_by_phoneme');

% Block boundary marker value inserted by the acquisition software
blockMarker = 800000;

% Which vowel is the standard / deviant in each block
stdPhonemeOfBlock = containers.Map({'A','B','C','D'}, {'bit','beat','bet','bat'});
devPhonemeOfBlock = containers.Map({'A','B','C','D'}, {'beat','bit','bat','bet'});
phonemes          = {'bit','beat','bet','bat'};

% Human-readable labels used in plot titles
phonemeLabels = containers.Map({'bit','beat','bet','bat'}, ...
                                {'/I/ (bit)', '/i:/ (beat)', '/e/ (bet)', '/ae/ (bat)'});

% Epoch window (seconds) and baseline (milliseconds)
epochWindow    = [-0.2  0.8];
baselineWindow = [-200  0];

% Frontocentral channels to include in the statistical test
chanLabels = {'F3','F4','FC3','FCz','FC4','C3','Cz','C4', 'F7', 'F8', 'Fp1', 'Fp2'};

% Analysis time window (seconds) for the permutation test
tWin = [0  0.8];

% Number of permutations (10000 recommended for publication)
nPerms = 10000;

% Cluster-forming and cluster-level significance threshold
clAlpha = 0.05;

% Colour scheme for grand-average plots (one colour per phoneme)
phonemeColours = containers.Map({'bit','beat','bet','bat'}, ...
                                 {[0 0 0.8], [0.2 0.4 1], [0.8 0 0], [1 0.4 0.2]});

%% -----------------------------------------------------------------------
%% Setup

if ~exist(outDir, 'dir'), mkdir(outDir); end

% Discover files and split into pre / post
allFiles = [];
for p = 1:numel(filePatterns)
    allFiles = [allFiles; dir(fullfile(baseDir, filePatterns{p}))]; %#ok<AGROW>
end
isPost    = contains({allFiles.name}, '_post');
preFiles  = allFiles(~isPost);
postFiles = allFiles(isPost);

getSubjID = @(n) regexp(n, '^([^_]+)', 'match', 'once');
preIDs    = cellfun(getSubjID, {preFiles.name},  'UniformOutput', false);
postIDs   = cellfun(getSubjID, {postFiles.name}, 'UniformOutput', false);

% Only subjects who completed BOTH sessions
completerIDs = intersect(preIDs, postIDs, 'stable');
fprintf('%d subjects completed both sessions — used for all phoneme tests.\n', numel(completerIDs));

preTbl  = readtable(preBlockTable,  'TextType', 'string');
postTbl = readtable(postBlockTable, 'TextType', 'string');

sessions = struct( ...
    'label', {'pre', 'post'}, ...
    'files', {preFiles, postFiles}, ...
    'tbl',   {preTbl, postTbl});

%% -----------------------------------------------------------------------
%% EXTRACTION: per-subject, per-session, per-phoneme standard/deviant ERPs

for sess = 1:numel(sessions)

    sessLabel = sessions(sess).label;
    pairTbl   = sessions(sess).tbl;
    files     = sessions(sess).files;
    fprintf('\n========== SESSION: %s ==========\n', upper(sessLabel));

    for f = 1:numel(files)

        fileName  = files(f).name;
        [~, baseName] = fileparts(fileName);
        subjectID = getSubjID(baseName);

        if ~ismember(subjectID, completerIDs)
            continue;   % skip subjects without data for both sessions
        end

        fprintf('\n  ----- %s -----\n', subjectID);

        rowIdx = find(pairTbl.Participant_ID == string(subjectID));
        if isempty(rowIdx)
            warning('No block order found for %s (%s) — skipping.', subjectID, sessLabel);
            continue;
        end
        blockLetters = [pairTbl.Block_1(rowIdx), pairTbl.Block_2(rowIdx), ...
                        pairTbl.Block_3(rowIdx), pairTbl.Block_4(rowIdx)];

        EEG = pop_loadset('filename', fileName, 'filepath', files(f).folder);
        EEG = eeg_checkset(EEG);
        EEG = pop_eegfiltnew(EEG, 'hicutoff', 30);

        allTypes = {EEG.event.type};
        isBlock  = cellfun(@(t) isequal(t, blockMarker) || ...
                                isequal(t, num2str(blockMarker)) || ...
                                (ischar(t) && strcmp(strtrim(t), num2str(blockMarker))), allTypes);
        blockLats = [EEG.event(isBlock).latency];

        if numel(blockLats) == 3
            blockBounds = [1, blockLats, EEG.pnts + 1];
        elseif numel(blockLats) == 4
            blockBounds = [blockLats, EEG.pnts + 1];
        else
            warning('%s (%s): expected 3 or 4 block markers, found %d — skipping.', ...
                    subjectID, sessLabel, numel(blockLats));
            continue;
        end

        % Per-phoneme role accumulators for this subject/session
        roleERP = struct();
        for ph = 1:numel(phonemes)
            roleERP.(phonemes{ph}).STD = [];
            roleERP.(phonemes{ph}).DEV = [];
        end
        times = []; chanlocs = [];

        for b = 1:4

            pb = char(blockLetters(b));
            if ~isKey(stdPhonemeOfBlock, pb)
                fprintf('    Block %d: unrecognised letter "%s" — skipping.\n', b, pb);
                continue;
            end

            EEG_blk = pop_select(EEG, 'point', [blockBounds(b), blockBounds(b+1)-1]);

            try
                EEG_s = pop_epoch(EEG_blk, {1}, epochWindow);
                EEG_s = pop_rmbase(EEG_s, baselineWindow);
            catch
                EEG_s = [];
            end
            try
                EEG_d = pop_epoch(EEG_blk, {2}, epochWindow);
                EEG_d = pop_rmbase(EEG_d, baselineWindow);
            catch
                EEG_d = [];
            end

            if isempty(EEG_s) || size(EEG_s.data,3) == 0 || ...
               isempty(EEG_d) || size(EEG_d.data,3) == 0
                fprintf('    Block %d (%s): no epochs found — skipping.\n', b, pb);
                continue;
            end

            times    = EEG_s.times;
            chanlocs = EEG_s.chanlocs;

            stdPh = stdPhonemeOfBlock(pb);
            devPh = devPhonemeOfBlock(pb);

            roleERP.(stdPh).STD = mean(EEG_s.data, 3);
            roleERP.(stdPh).nSTD = size(EEG_s.data, 3);
            roleERP.(devPh).DEV = mean(EEG_d.data, 3);
            roleERP.(devPh).nDEV = size(EEG_d.data, 3);

        end  % block loop

        if isempty(times), continue; end

        % Save one file per phoneme with both roles for this subject/session
        for ph = 1:numel(phonemes)
            phoneme = phonemes{ph};
            r = roleERP.(phoneme);
            if isempty(r.STD) || isempty(r.DEV)
                fprintf('    %s: missing standard or deviant role data — skipping.\n', phoneme);
                continue;
            end
            d.erpSTD   = r.STD;
            d.erpDEV   = r.DEV;
            d.times    = times;
            d.chanlocs = chanlocs;
            d.nSTD     = r.nSTD;
            d.nDEV     = r.nDEV;
            d.subject  = subjectID;
            d.session  = sessLabel;
            d.phoneme  = phoneme;
            save(fullfile(outDir, sprintf('%s_%s_role_%s.mat', subjectID, sessLabel, phoneme)), 'd');
            fprintf('    %-5s: %d STD-role trials | %d DEV-role trials\n', phoneme, r.nSTD, r.nDEV);
        end

    end  % subject loop
end  % session loop

%% -----------------------------------------------------------------------
%% Build channel neighbour structure from actual electrode positions

repList = dir(fullfile(outDir, sprintf('*_pre_role_%s.mat', phonemes{1})));
if isempty(repList)
    error('No extracted role files found in:\n  %s', outDir);
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
fprintf('\nNeighbours computed for %d channels.\n', numel(d_rep.chanlocs));

%% -----------------------------------------------------------------------
%% Main loop: one cluster permutation test per session x phoneme

results = struct();

for sess = 1:numel(sessions)

    sessLabel = sessions(sess).label;

    for ph = 1:numel(phonemes)

        phoneme    = phonemes{ph};
        phonLabel  = phonemeLabels(phoneme);
        fprintf('\n======= SESSION %s — PHONEME: %s (%s) =======\n', upper(sessLabel), phoneme, phonLabel);

        roleFiles = dir(fullfile(outDir, sprintf('*_%s_role_%s.mat', sessLabel, phoneme)));
        subjIDs   = cellfun(@(n) regexp(n, sprintf('^(.+?)_%s_role', sessLabel), 'tokens', 'once'), ...
                             {roleFiles.name}, 'UniformOutput', false);
        subjIDs   = [subjIDs{:}];

        nSubj = numel(subjIDs);
        fprintf('  %d subjects with both roles.\n', nSubj);
        if nSubj < 3
            warning('Too few subjects for %s / %s — skipping.', sessLabel, phoneme);
            continue;
        end

        ftSTD = {};
        ftDEV = {};
        for s = 1:nSubj
            tmp = load(fullfile(outDir, roleFiles(s).name), 'd');
            ftSTD{end+1} = toFieldTrip(tmp.d, 'STD'); %#ok<SAGROW>
            ftDEV{end+1} = toFieldTrip(tmp.d, 'DEV'); %#ok<SAGROW>
        end

        cfg                  = [];
        cfg.method           = 'montecarlo';
        cfg.statistic        = 'depsamplesT';
        cfg.correctm         = 'cluster';
        cfg.clusteralpha     = clAlpha;
        cfg.clusterstatistic = 'maxsum';
        cfg.tail             = 0;
        cfg.clustertail      = 0;
        cfg.correcttail      = 'prob';
        cfg.alpha            = clAlpha;
        cfg.numrandomization = nPerms;
        cfg.neighbours       = neighbours;
        cfg.channel          = {'all'};
        cfg.latency          = tWin;

        % Row 1 (ivar) — role: 1 = deviant, 2 = standard
        % Row 2 (uvar) — subject: 1..nSubj, repeated for each role
        cfg.design = [ones(1,nSubj), 2*ones(1,nSubj); ...
                      1:nSubj,       1:nSubj         ];
        cfg.ivar   = 1;
        cfg.uvar   = 2;

        fprintf('  Running permutation test (%d randomisations)...\n', nPerms);
        stat = ft_timelockstatistics(cfg, ftDEV{:}, ftSTD{:});

        printClusters(stat, phonLabel, sessLabel, nSubj, clAlpha);
        plotTStatistic(stat, phonLabel, sessLabel, nSubj, clAlpha);
        plotGrandAverage(ftSTD, ftDEV, stat, phoneme, phonLabel, sessLabel, nSubj, chanLabels, clAlpha, phonemeColours);

        results.(sessLabel).(phoneme)               = stat;
        results.(sessLabel).([phoneme '_subjects']) = subjIDs;

    end  % phoneme loop
end  % session loop

save(fullfile(outDir, 'ClusterPerm_StdVsDev_byPhoneme_FT.mat'), 'results');
fprintf('\nResults saved to:\n  %s\n', fullfile(outDir, 'ClusterPerm_StdVsDev_byPhoneme_FT.mat'));

%% -----------------------------------------------------------------------
%% Main loop 2: Pre vs Post comparison of the MMN difference wave (DEV-STD)
%% per phoneme

resultsPrePost = struct();

for ph = 1:numel(phonemes)

    phoneme   = phonemes{ph};
    phonLabel = phonemeLabels(phoneme);
    fprintf('\n======= PRE vs POST MMN — PHONEME: %s (%s) =======\n', phoneme, phonLabel);

    preRoleFiles  = dir(fullfile(outDir, sprintf('*_pre_role_%s.mat',  phoneme)));
    postRoleFiles = dir(fullfile(outDir, sprintf('*_post_role_%s.mat', phoneme)));

    preIDsPh  = cellfun(@(n) regexp(n, '^(.+?)_pre_role',  'tokens', 'once'), {preRoleFiles.name},  'UniformOutput', false);
    postIDsPh = cellfun(@(n) regexp(n, '^(.+?)_post_role', 'tokens', 'once'), {postRoleFiles.name}, 'UniformOutput', false);
    preIDsPh  = [preIDsPh{:}];
    postIDsPh = [postIDsPh{:}];

    pairedIDs = intersect(preIDsPh, postIDsPh, 'stable');
    nPaired   = numel(pairedIDs);
    fprintf('  %d subjects with both pre and post MMN data for %s.\n', nPaired, phoneme);

    if nPaired < 3
        warning('Too few paired subjects for phoneme %s — skipping.', phoneme);
        continue;
    end

    ftPre  = {};
    ftPost = {};
    for s = 1:nPaired
        fPre  = fullfile(outDir, sprintf('%s_pre_role_%s.mat',  pairedIDs{s}, phoneme));
        fPost = fullfile(outDir, sprintf('%s_post_role_%s.mat', pairedIDs{s}, phoneme));
        tmp = load(fPre,  'd');  ftPre{end+1}  = toFieldTripDiff(tmp.d);  %#ok<SAGROW>
        tmp = load(fPost, 'd');  ftPost{end+1} = toFieldTripDiff(tmp.d);  %#ok<SAGROW>
    end

    cfg                  = [];
    cfg.method           = 'montecarlo';
    cfg.statistic        = 'depsamplesT';
    cfg.correctm         = 'cluster';
    cfg.clusteralpha     = clAlpha;
    cfg.clusterstatistic = 'maxsum';
    cfg.tail             = 0;
    cfg.clustertail      = 0;
    cfg.correcttail      = 'prob';
    cfg.alpha            = clAlpha;
    cfg.numrandomization = nPerms;
    cfg.neighbours       = neighbours;
    cfg.channel          = {'all'};
    cfg.latency          = tWin;

    % Row 1 (ivar) — session: 1 = post, 2 = pre
    % Row 2 (uvar) — subject: 1..nPaired, repeated for each session
    cfg.design = [ones(1,nPaired), 2*ones(1,nPaired); ...
                  1:nPaired,       1:nPaired         ];
    cfg.ivar   = 1;
    cfg.uvar   = 2;

    fprintf('  Running permutation test (%d randomisations)...\n', nPerms);
    stat = ft_timelockstatistics(cfg, ftPost{:}, ftPre{:});

    printClustersPrePost(stat, phonLabel, nPaired, clAlpha);
    plotTStatisticPrePost(stat, phonLabel, nPaired, clAlpha);
    plotGrandAveragePrePost(ftPre, ftPost, stat, phoneme, phonLabel, nPaired, chanLabels, clAlpha, phonemeColours);

    resultsPrePost.(phoneme)               = stat;
    resultsPrePost.([phoneme '_subjects']) = pairedIDs;

end  % phoneme loop

save(fullfile(outDir, 'ClusterPerm_PreVsPost_byPhoneme_FT.mat'), 'resultsPrePost');
fprintf('\nResults saved to:\n  %s\n', fullfile(outDir, 'ClusterPerm_PreVsPost_byPhoneme_FT.mat'));

%% -----------------------------------------------------------------------
%% LOCAL FUNCTIONS
%% -----------------------------------------------------------------------

function tl = toFieldTrip(d, role)
% Convert an EEGLAB-style standard/deviant role struct to a FieldTrip
% timelock structure for the requested role ('STD' or 'DEV').
    if strcmp(role, 'STD')
        tl.avg = d.erpSTD;
    else
        tl.avg = d.erpDEV;
    end
    tl.time   = d.times / 1000;   % ms to seconds
    tl.label  = normaliseLabels({d.chanlocs.labels}');
    tl.dimord = 'chan_time';
end

% -------------------------------------------------------------------------
function tl = toFieldTripDiff(d)
% Convert a standard/deviant role struct to a FieldTrip timelock structure
% holding the MMN difference wave (DEV-role minus STD-role).
    tl.avg    = d.erpDEV - d.erpSTD;
    tl.time   = d.times / 1000;   % ms to seconds
    tl.label  = normaliseLabels({d.chanlocs.labels}');
    tl.dimord = 'chan_time';
end

% -------------------------------------------------------------------------
function labels = normaliseLabels(labels)
% Convert all-caps z-suffix labels to mixed-case (e.g. FCZ -> FCz, CZ -> Cz).
    for k = 1:numel(labels)
        lbl = labels{k};
        if numel(lbl) >= 2 && lbl(end) == 'Z'
            lbl(end)  = 'z';
            labels{k} = lbl;
        end
    end
end

% -------------------------------------------------------------------------
function printClusters(stat, phonLabel, sessLabel, nSubj, clAlpha)
% Print a summary of significant clusters to the console.
    fprintf('\n  ---- %s — Standard vs Deviant: %s  (N = %d) ----\n', upper(sessLabel), phonLabel, nSubj);
    fprintf('  Negative clusters (deviant-role MORE negative = genuine context effect):\n');
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
    fprintf('  Positive clusters (deviant-role LESS negative = unexpected direction):\n');
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
function plotTStatistic(stat, phonLabel, sessLabel, nSubj, clAlpha)
% Plot the t-statistic over time, averaged only across channels belonging
% to a significant cluster (falls back to all channels if none is
% significant), with significant cluster windows shaded. Green = negative
% cluster, orange = positive.
    sigChanIdx = significantClusterChannels(stat, clAlpha);
    if isempty(sigChanIdx)
        sigChanIdx = 1:size(stat.stat, 1);
    end
    tMean = mean(stat.stat(sigChanIdx, :), 1);
    tMs   = stat.time * 1000;

    try
        tCrit = tinv(1 - clAlpha/2, nSubj - 1);
    catch
        tCrit = 2;
        warning('tinv not available — using tCrit = 2 as fallback.');
    end
    yLim = max(max(abs(tMean)) * 1.4, tCrit * 1.4);

    figure('Color', 'white');
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
    xlabel('Time (ms)', 'FontName', 'Arial', 'FontSize', 10, 'FontWeight', 'bold');
    ylabel('t-statistic (channel average)', 'FontName', 'Arial', 'FontSize', 10, 'FontWeight', 'bold');
    set(gca, 'FontName', 'Arial', 'FontSize', 9);
    grid on; hold off;
end

% -------------------------------------------------------------------------
function plotGrandAverage(ftSTD, ftDEV, stat, phoneme, phonLabel, sessLabel, nSubj, chanLabels, clAlpha, colours)
% Grand-average ERP for standard-role vs deviant-role with cluster windows shaded.
    nCh = size(ftSTD{1}.avg, 1);
    nT  = size(ftSTD{1}.avg, 2);
    stackSTD = zeros(nCh, nT, nSubj);
    stackDEV = zeros(nCh, nT, nSubj);
    for s = 1:nSubj
        stackSTD(:,:,s) = ftSTD{s}.avg;
        stackDEV(:,:,s) = ftDEV{s}.avg;
    end
    gaSTD = mean(stackSTD, 3);
    gaDEV = mean(stackDEV, 3);

    times = ftSTD{1}.time * 1000;
    fcIdx = find(ismember(upper(ftSTD{1}.label), upper(chanLabels)));
    stdMean = mean(gaSTD(fcIdx, :), 1);
    devMean = mean(gaDEV(fcIdx, :), 1);

    baseCol = colours(phoneme);
    yLim    = max(max(abs([stdMean, devMean]))) * 1.4;
    if yLim == 0, yLim = 1; end

    figure('Color', 'white');
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

    hStd = plot(times, stdMean, '-',  'Color', baseCol, 'LineWidth', 2);
    hDev = plot(times, devMean, '--', 'Color', baseCol, 'LineWidth', 2);
    xline(0,'k:','LineWidth',1); yline(0,'k:','LineWidth',1);
    xlim([times(1) times(end)]); ylim([-yLim yLim]);
    xlabel('Time (ms)', 'FontName', 'Arial', 'FontSize', 10, 'FontWeight', 'bold');
    ylabel('Amplitude (\muV)', 'FontName', 'Arial', 'FontSize', 10, 'FontWeight', 'bold');
    set(gca, 'FontName', 'Arial', 'FontSize', 9);
    legend([hStd, hDev], {sprintf('Standard (n=%d)', nSubj), sprintf('Deviant (n=%d)', nSubj)}, 'Location','best','FontName', 'Arial', 'FontSize', 9, 'Box', 'off');
    grid off; hold off;
end

% -------------------------------------------------------------------------
function shadeBand(tMs, yLim, colour)
% Shade a time band with the given colour.
    if isempty(tMs), return; end
    patch([tMs(1) tMs(end) tMs(end) tMs(1)], [-yLim -yLim yLim yLim], ...
          colour, 'FaceAlpha', 0.20, 'EdgeColor', 'none');
end

% -------------------------------------------------------------------------
function chanIdx = significantClusterChannels(stat, clAlpha)
% Union of channel indices belonging to any negative/positive cluster that
% is significant at the cluster-level alpha threshold (empty if none is).
    chanIdx = [];
    if isfield(stat, 'negclusters')
        for k = 1:numel(stat.negclusters)
            if stat.negclusters(k).prob < clAlpha
                chanIdx = union(chanIdx, find(any(stat.negclusterslabelmat == k, 2)));
            end
        end
    end
    if isfield(stat, 'posclusters')
        for k = 1:numel(stat.posclusters)
            if stat.posclusters(k).prob < clAlpha
                chanIdx = union(chanIdx, find(any(stat.posclusterslabelmat == k, 2)));
            end
        end
    end
end

% -------------------------------------------------------------------------
function printClustersPrePost(stat, phonLabel, nPaired, clAlpha)
% Print a summary of significant clusters for the Pre vs Post MMN test.
    fprintf('\n  ---- Pre vs Post MMN: %s  (N = %d paired subjects) ----\n', phonLabel, nPaired);
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
function plotTStatisticPrePost(stat, phonLabel, nPaired, clAlpha)
% Plot the t-statistic over time for the Pre vs Post MMN test, averaged
% only across channels belonging to a significant cluster (falls back to
% all channels if none is significant), with significant cluster windows
% shaded. Green = negative (improvement), orange = positive (no
% improvement / regression).
    sigChanIdx = significantClusterChannels(stat, clAlpha);
    if isempty(sigChanIdx)
        sigChanIdx = 1:size(stat.stat, 1);
    end
    tMean = mean(stat.stat(sigChanIdx, :), 1);
    tMs   = stat.time * 1000;

    try
        tCrit = tinv(1 - clAlpha/2, nPaired - 1);
    catch
        tCrit = 2;
        warning('tinv not available — using tCrit = 2 as fallback.');
    end
    yLim = max(max(abs(tMean)) * 1.4, tCrit * 1.4);

    figure('Color', 'white');
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
    xlabel('Time (ms)', 'FontName', 'Arial', 'FontSize', 10, 'FontWeight', 'bold');
    ylabel('t-statistic (channel average)', 'FontName', 'Arial', 'FontSize', 10, 'FontWeight', 'bold');
    set(gca, 'FontName', 'Arial', 'FontSize', 9);
    grid on; hold off;
end

% -------------------------------------------------------------------------
function plotGrandAveragePrePost(ftPre, ftPost, stat, phoneme, phonLabel, nPaired, chanLabels, clAlpha, colours)
% Grand-average MMN difference wave for pre and post with cluster windows shaded.
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

    baseCol = colours(phoneme);
    yLim    = max(max(abs([preMean, postMean]))) * 1.4;
    if yLim == 0, yLim = 1; end

    figure('Color', 'white');
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
    xlabel('Time (ms)', 'FontName', 'Arial', 'FontSize', 10, 'FontWeight', 'bold');
    ylabel('Amplitude (\muV)', 'FontName', 'Arial', 'FontSize', 10, 'FontWeight', 'bold');
    set(gca, 'FontName', 'Arial', 'FontSize', 9);
    legend([hPre, hPost], {sprintf('Pre  (n=%d)', nPaired), sprintf('Post (n=%d)', nPaired)}, 'Location','best','FontName', 'Arial', 'FontSize', 9, 'Box', 'off');
    grid off; hold off;
end
