%% 04_MMN_by_contrast_pre_post.m
%
%  Fourth step in the MMN analysis pipeline.
%
%  Extracts subject-level and group-average MMN difference waves (DEV - STD)
%  separately for two phonetic contrasts, at pre-test and post-test.
%  Results are saved ready for the cluster-based permutation test in script 05.
%
%  The experiment used a block design in which each block contained stimuli
%  from one phonetic contrast. A CSV file (pairblock table) records which
%  contrast was presented in each block for each participant. This script
%  reads that table to route epochs to the correct contrast.
%
%  Contrasts (adapt labels and block letters in USER SETTINGS for your study)
%  --------------------------------------------------------------------------
%    Contrast 1 (AB) — block letters A or B  — e.g. /i:/ vs /I/
%    Contrast 2 (CD) — block letters C or D  — e.g. /ae/ vs /e/
%
%  Block boundary markers
%  ----------------------
%  Two conventions are handled automatically:
%    3 markers — each marker is a boundary between blocks
%    4 markers — each marker is the start of a block
%  If your data uses a different number of blocks or a different marker
%  value, adjust blockMarker, nBlocks, and the blockBounds logic below.
%
%  Outputs (saved to outDir)
%  -------------------------
%    <subjectID>_pre_MMN_AB.mat    subject-level MMN, pre-test,  contrast 1
%    <subjectID>_pre_MMN_CD.mat    subject-level MMN, pre-test,  contrast 2
%    <subjectID>_post_MMN_AB.mat   subject-level MMN, post-test, contrast 1
%    <subjectID>_post_MMN_CD.mat   subject-level MMN, post-test, contrast 2
%    GrandAverages_by_contrast.mat grand-average struct (pre/post x contrast)
%
%  Two figures are produced showing pre vs post grand-average MMN per contrast.

clear; eeglab;

%% -----------------------------------------------------------------------
%% USER SETTINGS — edit these before running

% Root directory of the project
baseDir = 'C:\MyProject';   % <-- change to your project folder

% Subfolders containing pre-test and post-test _GEDAIclean.set files
prePath  = fullfile(baseDir, 'Pre_test');
postPath = fullfile(baseDir, 'Post_test');

% CSV files mapping each participant to their block order
% Each file must have columns: Participant_ID, Block_1, Block_2, Block_3, Block_4
% Block values are single letters identifying the contrast (e.g. 'A','B','C','D')
preBlockTable  = fullfile(baseDir, 'pairblock_summary.csv');
postBlockTable = fullfile(baseDir, 'pairblock_summary_post.csv');

% Output directory for subject MMN files and grand averages
outDir = fullfile(baseDir, 'ERP_differences_by_contrast');

% File patterns to search for in each session folder
% Two patterns are listed here to handle variants in naming convention;
% add or remove patterns as needed for your dataset
filePatterns = {'*_GEDAIclean.set', '*_fixed_GEDAI_clean.set'};

% Block boundary marker value inserted by the acquisition software
blockMarker = 800000;

% Letters that identify each contrast in the pairblock table
contrastAB_letters = {'A', 'B'};   % contrast 1
contrastCD_letters = {'C', 'D'};   % contrast 2

% Epoch window (seconds) and baseline (milliseconds)
epochWindow    = [-0.2  0.8];
baselineWindow = [-200  0];

% Frontocentral channels used for the grand-average plots
fcChannels = {'F3','F4','FC3','FCZ','FC4','C3','CZ','C4'};

%% -----------------------------------------------------------------------
%% Setup

if ~exist(outDir, 'dir'), mkdir(outDir); end

preTbl  = readtable(preBlockTable,  'TextType', 'string');
postTbl = readtable(postBlockTable, 'TextType', 'string');

% Subject-level MMN containers
DIFF_pre_AB  = {};
DIFF_pre_CD  = {};
DIFF_post_AB = {};
DIFF_post_CD = {};

% Trial count log (printed as a summary table at the end)
trialLog = cell(0, 5);   % columns: subject, session, contrast, nSTD, nDEV

sessions = struct( ...
    'label', {'pre',    'post'   }, ...
    'path',  {prePath,  postPath }, ...
    'tbl',   {preTbl,   postTbl  });

%% -----------------------------------------------------------------------
%% Main loop: sessions x subjects

for sess = 1:numel(sessions)

    sessLabel = sessions(sess).label;
    filePath  = sessions(sess).path;
    pairTbl   = sessions(sess).tbl;

    % Collect files matching any of the specified patterns
    files = [];
    for p = 1:numel(filePatterns)
        files = [files; dir(fullfile(filePath, filePatterns{p}))]; %#ok<AGROW>
    end
    nFiles = numel(files);
    fprintf('\n========== SESSION: %s  (%d files) ==========\n', upper(sessLabel), nFiles);

    for f = 1:nFiles

        fileName  = files(f).name;
        [~, baseName] = fileparts(fileName);

        % Extract subject ID — everything up to the first underscore after
        % the participant prefix. Adjust the regex if your IDs differ.
        subjectID = regexp(baseName, '^([^_]+)', 'match', 'once');
        fprintf('\n  ----- %s -----\n', subjectID);

        %% Look up block order for this subject
        rowIdx = find(pairTbl.Participant_ID == string(subjectID));
        if isempty(rowIdx)
            warning('No block order found for %s (%s) — skipping.', subjectID, sessLabel);
            continue;
        end
        blockLetters = [pairTbl.Block_1(rowIdx), pairTbl.Block_2(rowIdx), ...
                        pairTbl.Block_3(rowIdx), pairTbl.Block_4(rowIdx)];

        %% Load and low-pass filter
        EEG = pop_loadset('filename', fileName, 'filepath', filePath);
        EEG = eeg_checkset(EEG);
        EEG = pop_eegfiltnew(EEG, 'hicutoff', 30);

        %% Locate block boundary markers and derive block boundaries
        %
        %  Two conventions are supported:
        %    3 markers — boundaries between 4 blocks; Block 1 starts at sample 1
        %    4 markers — start of each block; data before marker 1 is discarded
        %  If neither applies, the participant is skipped with a warning.

        allTypes  = {EEG.event.type};
        isBlock   = cellfun(@(t) isequal(t, blockMarker)  || ...
                                 isequal(t, num2str(blockMarker)) || ...
                                 (ischar(t) && strcmp(strtrim(t), num2str(blockMarker))), allTypes);
        blockLats = [EEG.event(isBlock).latency];

        if numel(blockLats) == 3
            % 3 boundary markers — Block 1 starts at sample 1
            blockBounds = [1, blockLats, EEG.pnts + 1];
        elseif numel(blockLats) == 4
            % 4 start-of-block markers — data before the first marker is discarded
            blockBounds = [blockLats, EEG.pnts + 1];
        else
            warning('%s (%s): expected 3 or 4 block markers, found %d — skipping.', ...
                    subjectID, sessLabel, numel(blockLats));
            continue;
        end

        %% Per-block epoch accumulators
        STD_AB = [];  DEV_AB = [];
        STD_CD = [];  DEV_CD = [];
        lastEEG_std = [];

        for b = 1:4

            pb = char(blockLetters(b));

            % Extract the continuous segment for this block
            EEG_blk = pop_select(EEG, 'point', [blockBounds(b), blockBounds(b+1)-1]);

            % Epoch standards and deviants; skip block if either is missing
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

            lastEEG_std = EEG_s;

            % Route epochs to the correct contrast bucket
            if ismember(pb, contrastAB_letters)
                STD_AB = cat(3, STD_AB, EEG_s.data);
                DEV_AB = cat(3, DEV_AB, EEG_d.data);
            elseif ismember(pb, contrastCD_letters)
                STD_CD = cat(3, STD_CD, EEG_s.data);
                DEV_CD = cat(3, DEV_CD, EEG_d.data);
            else
                fprintf('    Block %d: unrecognised letter "%s" — skipping.\n', b, pb);
            end

        end  % block loop

        if isempty(lastEEG_std), continue; end
        times    = lastEEG_std.times;
        chanlocs = lastEEG_std.chanlocs;

        %% Compute and save subject-level MMN (DEV - STD) per contrast

        if ~isempty(STD_AB) && ~isempty(DEV_AB)
            d.erp      = mean(DEV_AB, 3) - mean(STD_AB, 3);
            d.times    = times;
            d.chanlocs = chanlocs;
            d.nSTD     = size(STD_AB, 3);
            d.nDEV     = size(DEV_AB, 3);
            d.subject  = subjectID;
            if strcmp(sessLabel, 'pre')
                DIFF_pre_AB{end+1}  = d; %#ok<SAGROW>
            else
                DIFF_post_AB{end+1} = d; %#ok<SAGROW>
            end
            save(fullfile(outDir, sprintf('%s_%s_MMN_AB.mat', subjectID, sessLabel)), 'd');
            fprintf('    Contrast 1 (AB): %d STD  |  %d DEV trials\n', d.nSTD, d.nDEV);
            trialLog(end+1, :) = {subjectID, sessLabel, 'AB', d.nSTD, d.nDEV}; %#ok<SAGROW>
        end

        if ~isempty(STD_CD) && ~isempty(DEV_CD)
            d.erp      = mean(DEV_CD, 3) - mean(STD_CD, 3);
            d.times    = times;
            d.chanlocs = chanlocs;
            d.nSTD     = size(STD_CD, 3);
            d.nDEV     = size(DEV_CD, 3);
            d.subject  = subjectID;
            if strcmp(sessLabel, 'pre')
                DIFF_pre_CD{end+1}  = d; %#ok<SAGROW>
            else
                DIFF_post_CD{end+1} = d; %#ok<SAGROW>
            end
            save(fullfile(outDir, sprintf('%s_%s_MMN_CD.mat', subjectID, sessLabel)), 'd');
            fprintf('    Contrast 2 (CD): %d STD  |  %d DEV trials\n', d.nSTD, d.nDEV);
            trialLog(end+1, :) = {subjectID, sessLabel, 'CD', d.nSTD, d.nDEV}; %#ok<SAGROW>
        end

    end  % subject loop
end  % session loop

%% -----------------------------------------------------------------------
%% Trial count summary

fprintf('\n========== DEVIANT TRIAL COUNT SUMMARY ==========\n');
fprintf('%-12s  %-6s  %-10s  %6s  %6s\n', 'Subject','Session','Contrast','nSTD','nDEV');
fprintf('%s\n', repmat('-',1,48));
for r = 1:size(trialLog,1)
    fprintf('%-12s  %-6s  %-10s  %6d  %6d\n', trialLog{r,:});
end
fprintf('%s\n', repmat('-',1,48));

%% -----------------------------------------------------------------------
%% Grand averages across subjects for each condition

condNames = {'pre_AB',    'post_AB',    'pre_CD',    'post_CD'   };
condData  = {DIFF_pre_AB, DIFF_post_AB, DIFF_pre_CD, DIFF_post_CD};
GA = struct();

for c = 1:numel(condNames)

    cond     = condNames{c};
    subjData = condData{c};

    if isempty(subjData)
        warning('No data for condition %s.', cond);
        continue;
    end

    % Average over channels common to all subjects in this condition
    allLbls    = cellfun(@(x) {x.chanlocs.labels}, subjData, 'UniformOutput', false);
    commonLbls = allLbls{1};
    for s = 2:numel(allLbls)
        commonLbls = intersect(commonLbls, allLbls{s}, 'stable');
    end

    nS    = numel(subjData);
    nCh   = numel(commonLbls);
    nT    = numel(subjData{1}.times);
    stack = nan(nCh, nT, nS);

    for s = 1:nS
        lblsS    = {subjData{s}.chanlocs.labels};
        [~, idx] = ismember(commonLbls, lblsS);
        stack(:,:,s) = subjData{s}.erp(idx, :);
    end

    GA.(cond).erp      = mean(stack, 3, 'omitnan');
    GA.(cond).times    = subjData{1}.times;
    GA.(cond).labels   = commonLbls;
    [~, idx0]          = ismember(commonLbls, {subjData{1}.chanlocs.labels});
    GA.(cond).chanlocs = subjData{1}.chanlocs(idx0);
    GA.(cond).n        = nS;

end

save(fullfile(outDir, 'GrandAverages_by_contrast.mat'), 'GA');

%% -----------------------------------------------------------------------
%% Plots: pre vs post grand-average MMN per contrast

getFCmean = @(ga) mean( ...
    ga.erp(find(ismember(upper(ga.labels), upper(fcChannels))), :), 1);

% Figure 1 — Contrast 1 (AB)
figure; hold on;
if isfield(GA, 'pre_AB'),  plot(GA.pre_AB.times/1000,  getFCmean(GA.pre_AB),  'b-',  'LineWidth', 2); end
if isfield(GA, 'post_AB'), plot(GA.post_AB.times/1000, getFCmean(GA.post_AB), 'b--', 'LineWidth', 2); end
xline(0,'k:'); yline(0,'k:');
xlim([epochWindow(1) epochWindow(2)]); xlabel('Time (s)'); ylabel('Amplitude (\muV)');
legAB = {};
if isfield(GA,'pre_AB'),  legAB{end+1} = sprintf('Pre  (n = %d)', GA.pre_AB.n);  end
if isfield(GA,'post_AB'), legAB{end+1} = sprintf('Post (n = %d)', GA.post_AB.n); end
if ~isempty(legAB), legend(legAB); end
title('Grand-average MMN — Contrast 1 (AB): Pre vs Post');
grid on;

% Figure 2 — Contrast 2 (CD)
figure; hold on;
if isfield(GA, 'pre_CD'),  plot(GA.pre_CD.times/1000,  getFCmean(GA.pre_CD),  'r-',  'LineWidth', 2); end
if isfield(GA, 'post_CD'), plot(GA.post_CD.times/1000, getFCmean(GA.post_CD), 'r--', 'LineWidth', 2); end
xline(0,'k:'); yline(0,'k:');
xlim([epochWindow(1) epochWindow(2)]); xlabel('Time (s)'); ylabel('Amplitude (\muV)');
legCD = {};
if isfield(GA,'pre_CD'),  legCD{end+1} = sprintf('Pre  (n = %d)', GA.pre_CD.n);  end
if isfield(GA,'post_CD'), legCD{end+1} = sprintf('Post (n = %d)', GA.post_CD.n); end
if ~isempty(legCD), legend(legCD); end
title('Grand-average MMN — Contrast 2 (CD): Pre vs Post');
grid on;

fprintf('\nDone. Files saved to:\n  %s\n', outDir);
