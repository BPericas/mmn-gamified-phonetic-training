%% S06_MMN_pre_post_summary.m
%
%  Reports means and 95% CIs for the three cluster-permutation analyses
%  in this pipeline:
%
%    1) Overall MMN — DEV-STD pooled across all 30 recordings (15
%       subjects x pre/post), no pre/post distinction. Uses the
%       significant cluster from stat.mat (S04 Analysis 1, one-tailed).
%
%    2) Pre vs Post MMN — two-tailed, 15 paired subjects. Uses the
%       significant cluster from stats_prepostdiff.mat (S04 Analysis 2).
%       Reports Pre mean+CI, Post mean+CI, and the paired Post-Pre mean
%       change+CI (computed from each subject's own Post-Pre difference,
%       not from the two independent CIs).
%
%    3) Same-phoneme Standard vs Deviant (DEV_role - STD_role) — reported
%       separately for each phoneme x session, since that is how the
%       inference was actually run in S05 (8 separate cluster tests, not
%       one pooled test). Uses ClusterPerm_StdVsDev_byPhoneme_FT.mat.
%
%  Every mean/CI is computed over a FIXED set of channels and a fixed time
%  window (summaryChannels / summaryWinMs below), not from each test's own
%  significant cluster. The cluster-permutation stat files (stat.mat,
%  stats_prepostdiff.mat, ClusterPerm_StdVsDev_byPhoneme_FT.mat) are only
%  used to print whether each effect was statistically significant, as
%  context alongside the descriptive numbers.
%
%  Input:   ERP_differences/stat.mat, stats_prepostdiff.mat,
%           <subjectID>_ERPdiff.mat / <subjectID>_post_test_ERPdiff.mat
%           ERP_differences_by_phoneme/ClusterPerm_StdVsDev_byPhoneme_FT.mat,
%           <subjectID>_<session>_role_<phoneme>.mat
%  Output:  console summary
%           MMN_summary_all_analyses.csv — one row per reported value

clear;

%% -----------------------------------------------------------------------
%% USER SETTINGS — edit these before running

baseDir    = 'U:\SHaPS\Paul_Iverson\Bego_PhD\gamified_training';
erpDir     = fullfile(baseDir, 'ERP_differences');              % S04 output
phonemeDir = fullfile(baseDir, 'ERP_differences_by_phoneme');   % S05 output
clAlpha    = 0.05;   % cluster-level significance threshold used when the tests were run
phonemes   = {'bit','beat','bet','bat'};

% Descriptive MMN measurement settings
summaryChannels = {'F3','F4','FC3','FCz','FC4','C3','Cz','C4', ...
                   'F7','F8','Fp1','Fp2'};
summaryWinMs    = [144 232];
channelsStr     = strjoin(summaryChannels, ', ');

allRows = {};   % Analysis, Group, n, Mean, CI_low, CI_high, WindowStart_ms, WindowEnd_ms, Channels

%% =====================================================================
%% ANALYSIS 1: Overall MMN (no pre/post distinction)
%% =====================================================================
fprintf('\n============================================================\n');
fprintf('ANALYSIS 1: Overall MMN (DEV-STD), no pre/post distinction\n');
fprintf('============================================================\n');
try
    tmp = load(fullfile(erpDir, 'stat.mat'), 'stat');
    reportSignificance(tmp.stat, clAlpha, 'Overall MMN');
catch ME
    fprintf('(significance check skipped: %s)\n', ME.message);
end

try
    files = dir(fullfile(erpDir, 'PT*_ERPdiff.mat'));
    vals  = nan(numel(files), 1);
    for i = 1:numel(files)
        vals(i) = computeAmplitude(fullfile(files(i).folder, files(i).name), 'ERPdiff', summaryChannels, summaryWinMs);
    end

    [m, sd, ci] = meanCI95(vals);

    fprintf('DEV-STD (n=%d): Mean = %.3f uV, SD = %.3f, 95%% CI = [%.3f, %.3f]\n', numel(vals), m, sd, ci(1), ci(2));
    allRows(end+1,:) = {'1_Overall_MMN', 'DEV-STD', numel(vals), m, sd, ci(1), ci(2), summaryWinMs(1), summaryWinMs(2), channelsStr};
catch ME
    fprintf('Analysis 1 could not be completed: %s\n', ME.message);
end

%% =====================================================================
%% ANALYSIS 2: Pre vs Post MMN
%% =====================================================================
fprintf('\n============================================================\n');
fprintf('ANALYSIS 2: Pre vs Post MMN (two-tailed, paired subjects)\n');
fprintf('============================================================\n');
try
    statFile2 = fullfile(erpDir, 'stats_prepostdiff.mat');
    if ~isfile(statFile2)
        error('stats_prepostdiff.mat not found in:\n  %s\nRun S04 Analysis 2 first.', erpDir);
    end
    tmp = load(statFile2, 'stat2_prepostdiff');
    reportSignificance(tmp.stat2_prepostdiff, clAlpha, 'Pre vs Post MMN');
catch ME
    fprintf('(significance check skipped: %s)\n', ME.message);
end

try
    all_files  = dir(fullfile(erpDir, 'PT*_ERPdiff.mat'));
    post_files = all_files(contains({all_files.name}, '_post_test'));
    all_names  = {all_files.name};
    post_names = {post_files.name};
    pre_files  = all_files(~ismember(all_names, post_names));

    getID    = @(n) regexp(n, '^(PT[^_]+)', 'match', 'once');
    pre_ids  = cellfun(getID, {pre_files.name},  'UniformOutput', false);
    post_ids = cellfun(getID, {post_files.name}, 'UniformOutput', false);
    [matched_ids, pre_idx, post_idx] = intersect(pre_ids, post_ids, 'stable');
    nPaired = numel(matched_ids);
    if nPaired == 0
        error('No participants with both pre- and post-test ERPdiff files were found.');
    end

    preVals  = nan(nPaired, 1);
    postVals = nan(nPaired, 1);
    for i = 1:nPaired
        preFile  = pre_files(pre_idx(i));
        postFile = post_files(post_idx(i));
        preVals(i)  = computeAmplitude(fullfile(preFile.folder,  preFile.name),  'ERPdiff', summaryChannels, summaryWinMs);
        postVals(i) = computeAmplitude(fullfile(postFile.folder, postFile.name), 'ERPdiff', summaryChannels, summaryWinMs);
    end
    diffVals = postVals - preVals;   % per-subject paired difference

    [mPre,  sdPre,  ciPre]  = meanCI95(preVals);
    [mPost, sdPost, ciPost] = meanCI95(postVals);
    [mDiff, sdDiff, ciDiff] = meanCI95(diffVals);

    fprintf('Pre-test   (n=%d): Mean = %.3f uV, SD = %.3f, 95%% CI = [%.3f, %.3f]\n', ...
    nPaired, mPre, sdPre, ciPre(1), ciPre(2));

    fprintf('Post-test  (n=%d): Mean = %.3f uV, SD = %.3f, 95%% CI = [%.3f, %.3f]\n', ...
    nPaired, mPost, sdPost, ciPost(1), ciPost(2));

    fprintf('Post-Pre   (n=%d): Mean change = %.3f uV, SD = %.3f, 95%% CI = [%.3f, %.3f]\n', ...
    nPaired, mDiff, sdDiff, ciDiff(1), ciDiff(2));

    allRows(end+1,:) = {'2_Pre_vs_Post', 'Pre', nPaired, ...
    mPre, sdPre, ciPre(1), ciPre(2), ...
    summaryWinMs(1), summaryWinMs(2), channelsStr};

    allRows(end+1,:) = {'2_Pre_vs_Post', 'Post', nPaired, ...
    mPost, sdPost, ciPost(1), ciPost(2), ...
    summaryWinMs(1), summaryWinMs(2), channelsStr};

    allRows(end+1,:) = {'2_Pre_vs_Post', 'Post-Pre', nPaired, ...
    mDiff, sdDiff, ciDiff(1), ciDiff(2), ...
    summaryWinMs(1), summaryWinMs(2), channelsStr};
catch ME
    fprintf('Analysis 2 could not be completed: %s\n', ME.message);
end

%% =====================================================================
%% ANALYSIS 3: Same-phoneme Standard vs Deviant
%% =====================================================================
fprintf('\n============================================================\n');
fprintf('ANALYSIS 3: Same-phoneme Standard vs Deviant (DEV_role - STD_role)\n');
fprintf('============================================================\n');
try
    statFile3 = fullfile(phonemeDir, 'ClusterPerm_StdVsDev_byPhoneme_FT.mat');
    if ~isfile(statFile3)
        error('ClusterPerm_StdVsDev_byPhoneme_FT.mat not found in:\n  %s\nRun S05 first.', phonemeDir);
    end
    tmp     = load(statFile3, 'results');
    results = tmp.results;

    sessions = {'pre', 'post'};
    % Store subject-level values so that Post-Pre differences can be calculated within each phoneme after processing both sessions.
    phonemeVals = struct();
    for s = 1:numel(sessions)
        sessLabel = sessions{s};
        for p = 1:numel(phonemes)
            phoneme = phonemes{p};
            if ~isfield(results, sessLabel) || ~isfield(results.(sessLabel), phoneme)
                fprintf('%-5s / %-4s: no stat found — skipping.\n', phoneme, sessLabel);
                continue;
            end
            stat3   = results.(sessLabel).(phoneme);
            subjIDs = results.(sessLabel).([phoneme '_subjects']);
            nSubj   = numel(subjIDs);

            reportSignificance(stat3, clAlpha, sprintf('%s / %s', phoneme, sessLabel));

            vals = nan(nSubj, 1);
            for i = 1:nSubj
                roleFile = fullfile(phonemeDir, sprintf('%s_%s_role_%s.mat', subjIDs{i}, sessLabel, phoneme));
                vals(i)  = computeRoleDiffAmplitude(roleFile, summaryChannels, summaryWinMs);
            end

            [m, sd, ci] = meanCI95(vals);
            fprintf('%-5s / %-4s (n=%d): DEV-STD Mean = %.3f uV, SD = %.3f, 95%% CI = [%.3f, %.3f]\n', phoneme, sessLabel, nSubj, m, sd, ci(1), ci(2));
            allRows(end+1,:) = {'3_StdVsDev', sprintf('%s_%s', phoneme, sessLabel), nSubj, m, sd, ci(1), ci(2), summaryWinMs(1), summaryWinMs(2), channelsStr};

    % Store values and participant IDs for paired Pre vs Post calculation
    phonemeVals.(phoneme).(sessLabel).ids  = subjIDs;
    phonemeVals.(phoneme).(sessLabel).vals = vals;
        end
    end
    %% Paired Post-Pre differences within each phoneme
fprintf('\nPaired Post-Pre differences within phoneme:\n');

for p = 1:numel(phonemes)

    phoneme = phonemes{p};

    % Make sure both sessions are available
    if ~isfield(phonemeVals, phoneme) || ...
       ~isfield(phonemeVals.(phoneme), 'pre') || ...
       ~isfield(phonemeVals.(phoneme), 'post')

        fprintf('%-5s: Pre and/or Post data unavailable - skipping paired difference.\n', ...
            phoneme);
        continue;
    end

    preIDs   = phonemeVals.(phoneme).pre.ids;
    postIDs  = phonemeVals.(phoneme).post.ids;

    preValsP  = phonemeVals.(phoneme).pre.vals;
    postValsP = phonemeVals.(phoneme).post.vals;

    % Match participants across Pre and Post
    [matchedIDs, preIdx, postIdx] = intersect(preIDs, postIDs, 'stable');

    nPairedPhoneme = numel(matchedIDs);

    if nPairedPhoneme == 0
        fprintf('%-5s: No participants with both Pre and Post data.\n', ...
            phoneme);
        continue;
    end

    % Participant-level paired differences
    preMatched  = preValsP(preIdx);
    postMatched = postValsP(postIdx);
    diffValsP   = postMatched - preMatched;

    % Mean, SD and CI of paired change scores
    [mDiffP, sdDiffP, ciDiffP] = meanCI95(diffValsP);

    fprintf('%-5s Post-Pre (n=%d): Mean change = %.3f uV, SD = %.3f, 95%% CI = [%.3f, %.3f]\n', ...
        phoneme, nPairedPhoneme, mDiffP, sdDiffP, ...
        ciDiffP(1), ciDiffP(2));

    allRows(end+1,:) = {'3_StdVsDev', ...
        sprintf('%s_Post-Pre', phoneme), ...
        nPairedPhoneme, mDiffP, sdDiffP, ...
        ciDiffP(1), ciDiffP(2), ...
        summaryWinMs(1), summaryWinMs(2), channelsStr};

end
catch ME
    fprintf('Analysis 3 could not be completed: %s\n', ME.message);
end

%% =====================================================================
%% Save combined results table
%% =====================================================================
if ~isempty(allRows)
   T = cell2table(allRows, 'VariableNames', ...
    {'Analysis','Group','n','Mean','SD','CI_low','CI_high', ...
     'WindowStart_ms','WindowEnd_ms','Channels'});
    outFile = fullfile(baseDir, 'MMN_summary_all_analyses.csv');
    writetable(T, outFile);
    fprintf('\nCombined summary saved to:\n  %s\n', outFile);  
end

%% -----------------------------------------------------------------------
%% LOCAL FUNCTIONS
%% -----------------------------------------------------------------------

function reportSignificance(stat, clAlpha, label)
% Print whether the test has any cluster significant at clAlpha, and its
% p-value/direction/time window, for context alongside the fixed-window
% descriptive statistics. Does not affect how the means/CIs are computed.
    bestProb = Inf; bestType = ''; bestIdx = NaN;

    if isfield(stat, 'negclusters')
        for k = 1:numel(stat.negclusters)
            if stat.negclusters(k).prob < bestProb
                bestProb = stat.negclusters(k).prob; bestType = 'neg'; bestIdx = k;
            end
        end
    end
    if isfield(stat, 'posclusters')
        for k = 1:numel(stat.posclusters)
            if stat.posclusters(k).prob < bestProb
                bestProb = stat.posclusters(k).prob; bestType = 'pos'; bestIdx = k;
            end
        end
    end

    if isinf(bestProb)
        fprintf('[%s] No clusters found.\n', label);
        return;
    end
    if bestProb >= clAlpha
        fprintf('[%s] Not significant (best cluster p = %.4f).\n', label, bestProb);
        return;
    end

    if strcmp(bestType, 'neg')
        labelmat = stat.negclusterslabelmat;
    else
        labelmat = stat.posclusterslabelmat;
    end
    timeIdx = any((labelmat == bestIdx) & logical(stat.mask), 1);
    winMs   = [stat.time(find(timeIdx,1,'first')), stat.time(find(timeIdx,1,'last'))] * 1000;
    fprintf('[%s] Significant %s cluster: p = %.4f, %.0f-%.0f ms.\n', label, bestType, bestProb, winMs(1), winMs(2));
end

function val = computeAmplitude(filePath, varName, sigChannels, timeWindowMs)
    S          = load(filePath, varName);
    erpStruct  = S.(varName);
    chanLabels = {erpStruct.chanlocs.labels};
    chanIdx    = ismember(chanLabels, sigChannels);
    timeIdx    = erpStruct.times >= timeWindowMs(1) & erpStruct.times <= timeWindowMs(2);
    val        = mean(mean(erpStruct.erp(chanIdx, timeIdx), 2), 1);
end

function val = computeRoleDiffAmplitude(filePath, sigChannels, timeWindowMs)
% MMN-style difference (DEV-role minus STD-role) for one subject/phoneme/session.
    S          = load(filePath, 'd');
    d          = S.d;
    chanLabels = {d.chanlocs.labels};
    chanIdx    = ismember(chanLabels, sigChannels);
    timeIdx    = d.times >= timeWindowMs(1) & d.times <= timeWindowMs(2);
    diffErp    = d.erpDEV - d.erpSTD;
    val        = mean(mean(diffErp(chanIdx, timeIdx), 2), 1);
end

function [m, sd, ci] = meanCI95(x)
x = x(~isnan(x));
n = numel(x);
m = mean(x);
sd = std(x);
sem = sd / sqrt(n);
tcrit = tinv(0.975, n - 1);
ci = [m - tcrit * sem, m + tcrit * sem];
end