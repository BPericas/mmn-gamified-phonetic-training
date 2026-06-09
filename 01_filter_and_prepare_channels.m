%% 01_filter_and_prepare_channels.m
%
%  First preprocessing step in the MMN analysis pipeline.
%
%  For each participant this script:
%    1. Loads the raw EEG file (.cdt from CURRY 8 or .set)
%    2. Applies a 0.1 Hz high-pass filter separately to each recording
%       block (blocks are separated by event marker 800000)
%    3. Saves the VEOG channel for later use as an artefact criterion
%    4. Removes the VEOG and trigger channels
%    5. Assigns standard 10-20 electrode coordinates
%    6. Standardises channel labels to mixed-case convention (e.g. Fz, Cz)
%    7. Saves the processed data as a .set file ready for GEDAI
%
%  Input:   raw .cdt (or .set) files in the folder specified by dataDir
%  Output:  <subjectID>_filtnotrig.set  — filtered, channels trimmed
%           <subjectID>_VEOG.mat        — VEOG time course saved separately
%           <subjectID>_filtnotrig.mat  — MATLAB snapshot of processed EEG
%
%  Notes
%  -----
%  - Filtering is applied per block rather than across the full continuous
%    recording to avoid filter edge artefacts at block boundaries.
%  - The VEOG channel is saved before removal so it can be used post-GEDAI
%    to identify eye-blink components.
%  - The trigger channel carries no neural signal and is discarded.
%  - Channel labels are converted from all-caps (NeuroScan convention) to
%    mixed-case (GEDAI/FieldTrip convention), e.g. FZ -> Fz, FCZ -> FCz.
%  - Files that have already been processed are skipped automatically.

clear;
eeglab;

%% -----------------------------------------------------------------------
%% USER SETTINGS — edit these before running

% Folder containing the raw EEG files
dataDir = pwd;   % defaults to current folder; replace with an explicit path if needed
% Example: dataDir = 'C:\MyProject\RawData';

% File pattern to match your participant files
% Use * as a wildcard. Examples:
%   'sub-*.cdt'      for BIDS-style naming
%   'Participant*.cdt'
%   '*.cdt'          to process all .cdt files in the folder
filePattern = '*.cdt';

% Index of the VEOG channel (used if the label 'VEOG' is not found)
veogChannelIndex = 31;

% Index of the trigger channel to remove (carries no neural signal)
triggerChannelIndex = 33;

%% -----------------------------------------------------------------------
%% Find participant files

files  = dir(fullfile(dataDir, filePattern));
nFiles = length(files);
fprintf('Found %d file(s) matching "%s".\n', nFiles, filePattern);
if nFiles == 0
    error('No files found. Check dataDir and filePattern in the USER SETTINGS section.');
end

%% -----------------------------------------------------------------------
%% Main loop — one iteration per participant

for f = 1:nFiles

    fileName  = files(f).name;
    [~, subjectID, ~] = fileparts(fileName);

    %% Skip already-processed files
    outMatFile = fullfile(files(f).folder, [subjectID '_filtnotrig.mat']);
    if exist(outMatFile, 'file')
        fprintf('\n>>> Skipping %s — already processed.\n', subjectID);
        continue;
    end

    fprintf('\n===== Processing %s (%d/%d) =====\n', subjectID, f, nFiles);

    %% STEP 1: Load raw EEG data
    [~, ~, ext] = fileparts(fileName);
    fprintf('  Loading %s ...\n', fileName);

    switch lower(ext)
        case '.cdt'
            % CURRY 8 format — try specifying version first, fall back if needed
            try
                EEG = pop_loadcurry(fullfile(files(f).folder, fileName), [], 'CurryVersion', 8);
            catch
                EEG = pop_loadcurry(fullfile(files(f).folder, fileName));
            end
        case '.set'
            EEG = pop_loadset('filename', fileName, 'filepath', files(f).folder);
        otherwise
            error('Unsupported file format: %s', ext);
    end

    EEG = eeg_checkset(EEG);
    fprintf('  Loaded: %d samples x %d channels (%.1f min @ %d Hz)\n', ...
        EEG.pnts, EEG.nbchan, EEG.pnts / EEG.srate / 60, EEG.srate);

    %% STEP 2: High-pass filter (0.1 Hz) applied per recording block
    %
    %  Recording blocks are separated by event marker 800000 inserted by
    %  CURRY 8. Filtering each block independently avoids edge artefacts
    %  that would arise from filtering across block boundaries.
    %  If your data does not use block markers, remove this block-splitting
    %  logic and apply pop_eegfiltnew directly to the full EEG structure.

    fprintf('  Identifying recording blocks (marker 800000)...\n');
    eventTypes  = [EEG.event.type];
    pauseEvents = find(eventTypes == 800000);

    blockStarts = 1;
    blockEnds   = [];
    for ei = 1:numel(pauseEvents)
        eventSample = round(EEG.event(pauseEvents(ei)).latency);
        blockEnds   = [blockEnds;   eventSample    ]; %#ok<AGROW>
        blockStarts = [blockStarts; eventSample + 1]; %#ok<AGROW>
    end
    blockEnds   = [blockEnds;   EEG.pnts];
    blockStarts = blockStarts(1:numel(blockEnds));

    fprintf('  Found %d block(s).\n', numel(blockStarts));

    for blk = 1:numel(blockStarts)
        s1 = blockStarts(blk);
        s2 = blockEnds(blk);
        fprintf('    Block %d: samples %d-%d (%.1f s)\n', blk, s1, s2, (s2-s1)/EEG.srate);

        EEG_blk      = EEG;
        EEG_blk.data = EEG.data(:, s1:s2);
        EEG_blk.pnts = s2 - s1 + 1;
        EEG_blk      = pop_eegfiltnew(EEG_blk, 'locutoff', 0.1);

        EEG.data(:, s1:s2) = EEG_blk.data;
    end

    EEG = eeg_checkset(EEG);

    %% STEP 3: Save VEOG channel before removal
    %
    %  The VEOG channel is saved here so it can be used after GEDAI to
    %  check whether any retained components correlate with eye movements.

    veog_idx = find(strcmpi({EEG.chanlocs.labels}, 'VEOG'));
    if isempty(veog_idx)
        veog_idx = veogChannelIndex;
        fprintf('  Warning: VEOG label not found — using channel index %d.\n', veogChannelIndex);
    end
    VEOG = double(EEG.data(veog_idx, :));
    save(fullfile(files(f).folder, [subjectID '_VEOG.mat']), 'VEOG', '-v7.3');
    fprintf('  VEOG saved (%d samples).\n', length(VEOG));

    %% STEP 4: Remove VEOG and trigger channels
    %
    %  Adjust veogChannelIndex and triggerChannelIndex in USER SETTINGS
    %  at the top of this script if your channel layout differs.

    EEG = pop_select(EEG, 'nochannel', [veogChannelIndex triggerChannelIndex]);
    EEG = eeg_checkset(EEG);

    %% STEP 5: Assign standard 10-20 electrode coordinates

    fprintf('  Assigning electrode coordinates...\n');
    EEG = pop_chanedit(EEG, 'lookup', ...
        fullfile(fileparts(which('eeglab.m')), ...
        'plugins', 'dipfit', 'standard_BESA', 'standard-10-5-cap385.elp'));
    EEG = eeg_checkset(EEG);

    %% STEP 6: Standardise channel labels to mixed-case convention
    %
    %  NeuroScan exports all-caps labels (FZ, FCZ, CZ). GEDAI and FieldTrip
    %  expect mixed-case labels (Fz, FCz, Cz). Remove this step if your
    %  system already exports mixed-case labels.

    for ch = 1:length(EEG.chanlocs)
        lbl = EEG.chanlocs(ch).labels;
        lbl = regexprep(lbl, '^FP', 'Fp');   % FP1 -> Fp1, FP2 -> Fp2
        lbl = regexprep(lbl, 'Z$',  'z');    % FZ  -> Fz,  CZ  -> Cz, FCZ -> FCz
        EEG.chanlocs(ch).labels = lbl;
    end
    fprintf('  Channel labels after standardisation:\n');
    disp({EEG.chanlocs.labels});

    %% STEP 7: Save processed data

    EEG = pop_saveset(EEG, ...
        'filename', [subjectID '_filtnotrig.set'], ...
        'filepath', files(f).folder);

    save(outMatFile, 'EEG', 'subjectID', '-v7.3');
    fprintf('  Saved: %s_filtnotrig.set\n', subjectID);

end

fprintf('\nDone. Processed %d file(s).\n', nFiles);
