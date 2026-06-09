%% 03_extract_ERPs_and_MMN.m
%
%  Third step in the MMN analysis pipeline.
%
%  For each participant this script:
%    1. Loads the artefact-cleaned EEG file produced by script 02
%    2. Applies a 30 Hz low-pass filter (standard for MMN ERP visualisation)
%    3. Epochs the data around standard (marker 1) and deviant (marker 2)
%       stimuli with a 200 ms pre-stimulus baseline
%    4. Applies baseline correction (-200 to 0 ms)
%    5. Computes subject-level ERPs by averaging across trials
%    6. Computes the MMN difference wave (DEV - STD) per subject and saves it
%    7. Computes grand-average ERPs and MMN across all subjects
%    8. Plots the group-average standard vs deviant ERP and the MMN
%
%  Input:   <subjectID>_GEDAIclean.set   (output of 02_run_GEDAI_artefact_removal.m)
%  Output:  ERP_differences/<subjectID>_ERPdiff.mat   — subject-level MMN
%                                                        difference wave
%
%  Trigger codes
%  -------------
%    1  = standard stimulus
%    2  = deviant stimulus
%  Change stdMarker and devMarker in USER SETTINGS if your codes differ.
%
%  Notes
%  -----
%  - Grand averages are computed over channels common to all subjects,
%    to handle cases where a channel was removed for one participant.
%  - The MMN difference wave is computed as DEV - STD. The MMN is a
%    negative deflection, so the difference wave should go negative
%    in the 100-250 ms window at frontocentral channels.

clear;
eeglab;

%% -----------------------------------------------------------------------
%% USER SETTINGS — edit these before running

% Folder containing the _GEDAIclean.set files
dataDir = pwd;   % defaults to current folder; replace with an explicit path if needed
% Example: dataDir = 'C:\MyProject\PreprocessedData';

% File pattern to match cleaned files from script 02
filePattern = '*_GEDAIclean.set';

% Output folder for subject-level difference waves
outDir = fullfile(dataDir, 'ERP_differences');

% Stimulus trigger codes
stdMarker = {1};   % standard stimulus
devMarker = {2};   % deviant stimulus

% Epoch window (seconds) and baseline period (milliseconds)
epochWindow   = [-0.2  0.8];
baselineWindow = [-200  0];

% Frontocentral channels used for grand-average plots
% Adjust to match your electrode layout
fcChannels = {'F3','F4','FC3','FCZ','FC4','C3','CZ','C4'};

%% -----------------------------------------------------------------------
%% Setup

if ~exist(outDir, 'dir')
    mkdir(outDir);
end

files  = [dir(fullfile(dataDir, filePattern))];
nFiles = numel(files);

if nFiles == 0
    error('No files matching "%s" found in:\n  %s\nCheck dataDir and filePattern.', filePattern, dataDir);
end

fprintf('Found %d file(s) to process.\n', nFiles);

% Containers for grand-average computation
allSTD = cell(nFiles, 1);
allDEV = cell(nFiles, 1);

%% -----------------------------------------------------------------------
%% Subject loop

for f = 1:nFiles

    fileName  = files(f).name;
    [~, baseName] = fileparts(fileName);
    subjectID = regexprep(baseName, '_GEDAIclean$', '');

    fprintf('\n===== Processing %s (%d/%d) =====\n', subjectID, f, nFiles);

    %% STEP 1: Load artefact-cleaned EEG
    EEG = pop_loadset('filename', fileName, 'filepath', files(f).folder);
    EEG = eeg_checkset(EEG);

    %% STEP 2: Low-pass filter at 30 Hz
    % Applied here rather than in script 02 so that GEDAI operates on
    % broadband data. 30 Hz is standard for MMN ERP visualisation.
    EEG = pop_eegfiltnew(EEG, 'hicutoff', 30);
    EEG = eeg_checkset(EEG);

    %% STEP 3: Epoch around standard and deviant stimuli
    EEG_std = pop_epoch(EEG, stdMarker, epochWindow);
    EEG_dev = pop_epoch(EEG, devMarker, epochWindow);

    %% STEP 4: Baseline correction
    EEG_std = pop_rmbase(EEG_std, baselineWindow);
    EEG_dev = pop_rmbase(EEG_dev, baselineWindow);

    % Verify that epochs were found
    assert(size(EEG_std.data, 3) > 0, 'No standard trials found for %s', subjectID);
    assert(size(EEG_dev.data, 3) > 0, 'No deviant trials found for %s',  subjectID);
    fprintf('  Trials: %d standard  |  %d deviant\n', ...
        size(EEG_std.data,3), size(EEG_dev.data,3));

    %% STEP 5: Compute subject-level ERPs (trial average)
    ERP_std = mean(EEG_std.data, 3);   % channels x time
    ERP_dev = mean(EEG_dev.data, 3);

    % Store for grand-average computation
    allSTD{f}.erp      = ERP_std;
    allSTD{f}.times    = EEG_std.times;
    allSTD{f}.chanlocs = EEG_std.chanlocs;
    allSTD{f}.ntrials  = size(EEG_std.data, 3);

    allDEV{f}.erp      = ERP_dev;
    allDEV{f}.times    = EEG_dev.times;
    allDEV{f}.chanlocs = EEG_dev.chanlocs;
    allDEV{f}.ntrials  = size(EEG_dev.data, 3);

    %% STEP 6: Compute and save subject-level MMN difference wave (DEV - STD)
    assert(isequal(size(ERP_dev), size(ERP_std)), 'Channel/time mismatch between STD and DEV ERPs');
    assert(isequal(EEG_dev.times, EEG_std.times), 'Time vectors differ between STD and DEV');

    ERPdiff.erp      = ERP_dev - ERP_std;
    ERPdiff.times    = EEG_std.times;
    ERPdiff.chanlocs = EEG_std.chanlocs;

    save(fullfile(outDir, [subjectID '_ERPdiff.mat']), 'ERPdiff');
    fprintf('  MMN difference wave saved.\n');

end

%% -----------------------------------------------------------------------
%% STEP 7: Grand-average ERP and MMN across all subjects
%
%  Averages are computed over channels common to all subjects, to handle
%  cases where a channel was removed for a specific participant.

fprintf('\nComputing grand averages...\n');

% Find common channel labels across all subjects
allLabels  = cellfun(@(x) {x.chanlocs.labels}, allSTD, 'UniformOutput', false);
commonLabels = allLabels{1};
for s = 2:numel(allLabels)
    commonLabels = intersect(commonLabels, allLabels{s}, 'stable');
end

nSubj = numel(allSTD);
nChan = numel(commonLabels);
nTime = numel(allSTD{1}.times);

STD_stack = nan(nChan, nTime, nSubj);
DEV_stack = nan(nChan, nTime, nSubj);

for s = 1:nSubj
    labels_s = {allSTD{s}.chanlocs.labels};
    [~, idx] = ismember(commonLabels, labels_s);
    STD_stack(:,:,s) = allSTD{s}.erp(idx, :);
    DEV_stack(:,:,s) = allDEV{s}.erp(idx, :);
end

GA_STD.erp = mean(STD_stack, 3, 'omitnan');
GA_DEV.erp = mean(DEV_stack, 3, 'omitnan');
GA_STD.times = allSTD{1}.times;
GA_DEV.times = allDEV{1}.times;
[~, idx0] = ismember(commonLabels, {allSTD{1}.chanlocs.labels});
GA_STD.chanlocs = allSTD{1}.chanlocs(idx0);
GA_DEV.chanlocs = GA_STD.chanlocs;

GA_DIFF.erp      = GA_DEV.erp - GA_STD.erp;
GA_DIFF.times    = GA_STD.times;
GA_DIFF.chanlocs = GA_STD.chanlocs;

%% -----------------------------------------------------------------------
%% STEP 8: Plot grand-average results

% Find frontocentral channel indices for plotting
gaLabels = {GA_STD.chanlocs.labels};
chanIdx  = find(ismember(upper(gaLabels), upper(fcChannels)));
if isempty(chanIdx)
    warning('None of the specified fcChannels found — plotting all channels.');
    chanIdx = 1:size(GA_STD.erp, 1);
end

STD_fc  = mean(GA_STD.erp(chanIdx,  :), 1);
DEV_fc  = mean(GA_DEV.erp(chanIdx,  :), 1);
DIFF_fc = mean(GA_DIFF.erp(chanIdx, :), 1);
times_s = GA_STD.times / 1000;   % ms to seconds

% Figure 1: Standard vs Deviant ERP
figure; hold on;
plot(times_s, STD_fc,  'LineWidth', 2);
plot(times_s, DEV_fc,  'LineWidth', 2);
xline(0, 'k:'); yline(0, 'k:');
xlim([epochWindow(1) epochWindow(2)]);
xlabel('Time (s)');
ylabel('Amplitude (\muV)');
legend({'Standard', 'Deviant'});
title(sprintf('Group-average ERP — frontocentral (N = %d)', nSubj));
grid on;

% Figure 2: MMN difference wave (DEV - STD)
figure; hold on;
plot(times_s, DIFF_fc, 'k', 'LineWidth', 2);
xline(0, 'k:'); yline(0, 'k:');
xlim([epochWindow(1) epochWindow(2)]);
xlabel('Time (s)');
ylabel('Amplitude (\muV)');
title(sprintf('Group-average MMN (DEV \x2212 STD) — frontocentral (N = %d)', nSubj));
grid on;

fprintf('\nDone. Subject MMN files saved to:\n  %s\n', outDir);
