%% 02_run_GEDAI_artefact_removal.m
%
%  Second preprocessing step in the MMN analysis pipeline.
%
%  Runs GEDAI (Generalised Eigenvalue Decomposition for Artefact
%  Identification) on all filtered files produced by script 01. After
%  GEDAI, a robust average rereference is applied using NoiseTools (nt_rereference).
%
%  Input:   <subjectID>_filtnotrig.set   (output of 01_filter_and_prepare_channels.m)
%           <subjectID>_VEOG.mat         (VEOG time course saved by script 01)
%  Output:  <subjectID>_GEDAIclean.set   — artefact-cleaned, rereferenced EEG
%           <subjectID>_GEDAIclean.mat   — MATLAB snapshot including GEDAI
%                                          quality metrics and VEOG correlations
%
%  Dependencies
%  ------------
%  - GEDAI toolbox  (https://github.com/NeuroEngUAB/GEDAI)
%  - NoiseTools     (http://audition.ens.fr/adc/NoiseTools/)
%
%  Notes
%  -----
%  - VEOG correlations are computed after GEDAI to quantify any residual
%    eye-blink contamination. These are saved for inspection but no
%    automatic exclusion is applied here.
%  - Files that have already been processed are skipped automatically.
%  - Parallel processing is disabled by default. If you have the Parallel
%    Computing Toolbox, set use_parallel = true for faster processing.

clear;
eeglab;

%% -----------------------------------------------------------------------
%% USER SETTINGS — edit these before running

% Folder containing the _filtnotrig.set files
dataDir = pwd;   % defaults to current folder; replace with an explicit path if needed
% Example: dataDir = 'C:\MyProject\PreprocessedData';

% File pattern to match processed files from script 01
% Change the prefix to match your participant naming convention
filePattern = '*_filtnotrig.set';

%% -----------------------------------------------------------------------
%% GEDAI parameters — adjust if needed

artifact_threshold  = 'auto';        % denoising strength: 'auto', 'auto+', 'auto-'
epoch_size_cycles   = 12;            % epoch size in wavelet cycles per frequency band
lowcut_frequency    = 0.5;           % exclude wavelet bands below this frequency (Hz)
ref_matrix_type     = 'interpolated';% leadfield type: 'interpolated' for non-standard locations
use_parallel        = false;         % set true if Parallel Computing Toolbox is available
visualize_artifacts = false;         % set true to inspect artefact components interactively
visualize_manifold  = false;         % set true to display the SENSAI manifold plot
ENOVA_threshold     = [];            % [] = no epoch rejection (epoching is done in script 03)

%% -----------------------------------------------------------------------
%% Find input files

files  = dir(fullfile(dataDir, filePattern));
nFiles = length(files);

if nFiles == 0
    error('No files matching "%s" found in:\n  %s\nCheck dataDir and filePattern.', filePattern, dataDir);
end

fprintf('Found %d file(s) to process.\n\n', nFiles);

%% -----------------------------------------------------------------------
%% Main loop — one iteration per participant

for f = 1:nFiles

    fileName = files(f).name;
    fileDir  = files(f).folder;

    % Extract subject ID by stripping the '_filtnotrig' suffix
    subjectID = regexp(fileName, '^(.+?)_filtnotrig', 'tokens', 'once');
    if isempty(subjectID)
        warning('Could not parse subject ID from filename: %s — skipping.', fileName);
        continue;
    end
    subjectID = subjectID{1};

    %% Skip already-processed files
    outMatFile = fullfile(fileDir, [subjectID '_GEDAIclean.mat']);
    if exist(outMatFile, 'file')
        fprintf('>>> Skipping %s — output already exists.\n', subjectID);
        continue;
    end

    fprintf('===== Processing %s (%d/%d) =====\n', subjectID, f, nFiles);

    %% STEP 1: Load filtered EEG data
    fprintf('  Loading %s ...\n', fileName);
    EEG = pop_loadset('filename', fileName, 'filepath', fileDir);
    EEG = eeg_checkset(EEG);
    fprintf('  Data: %d channels x %d samples (%.1f min @ %d Hz)\n', ...
        EEG.nbchan, EEG.pnts, EEG.pnts / EEG.srate / 60, EEG.srate);

    %% STEP 2: Load VEOG time course (saved by script 01)
    veogFile = fullfile(fileDir, [subjectID '_VEOG.mat']);
    if exist(veogFile, 'file')
        tmp     = load(veogFile, 'VEOG');
        VEOG    = tmp.VEOG(:);
        hasVEOG = true;
        fprintf('  VEOG loaded (%d samples).\n', length(VEOG));
        if length(VEOG) ~= EEG.pnts
            warning('VEOG length (%d) does not match EEG length (%d) — skipping VEOG analysis.', ...
                length(VEOG), EEG.pnts);
            hasVEOG = false;
        end
    else
        VEOG    = [];
        hasVEOG = false;
        fprintf('  Warning: no VEOG file found for %s.\n', subjectID);
    end

    %% STEP 3: Run GEDAI artefact removal
    fprintf('  Running GEDAI...\n');
    [EEG, ~, SENSAI_score, ~, ~, mean_ENOVA, ~, ~] = GEDAI( ...
        EEG, ...
        artifact_threshold, ...
        epoch_size_cycles, ...
        lowcut_frequency, ...
        ref_matrix_type, ...
        use_parallel, ...
        visualize_artifacts, ...
        ENOVA_threshold, ...
        [], ...              % signal_type: [] defaults to EEG
        visualize_manifold);

    fprintf('  SENSAI score: %.1f%%  |  Mean ENOVA: %.3f\n', SENSAI_score, mean_ENOVA);
    EEG = eeg_checkset(EEG);

    %% STEP 4: Post-GEDAI VEOG residual correlation check
    %
    %  Correlates the VEOG time course with each cleaned EEG channel.
    %  High correlations suggest residual eye-blink contamination.
    %  Results are saved for inspection — no automatic exclusion is applied.
    if hasVEOG
        veog_z     = zscore(VEOG);
        clean_data = double(EEG.data);
        veog_corr  = zeros(EEG.nbchan, 1);
        for ch = 1:EEG.nbchan
            r = corrcoef(clean_data(ch,:)', veog_z);
            veog_corr(ch) = r(1,2);
        end
        [max_r, max_ch] = max(abs(veog_corr));
        fprintf('  Max |VEOG corr| post-GEDAI: %.3f (channel %s)\n', ...
            max_r, EEG.chanlocs(max_ch).labels);
    else
        veog_corr = [];
    end

    %% STEP 5: Robust average rereference (NoiseTools)
    %
    %  nt_rereference computes a weighted average reference that is robust
    %  to channels with high noise, unlike a simple average reference.
    fprintf('  Applying robust average reference (NoiseTools)...\n');
    x        = double(EEG.data)';   % NoiseTools expects [samples x channels]
    x        = nt_rereference(x);
    EEG.data = x';                  % back to [channels x samples]
    EEG      = eeg_checkset(EEG);

    %% STEP 6: Save cleaned data
    fprintf('  Saving %s_GEDAIclean ...\n', subjectID);
    EEG = pop_saveset(EEG, ...
        'filename', [subjectID '_GEDAIclean.set'], ...
        'filepath', fileDir);

    save(outMatFile, 'EEG', 'subjectID', 'SENSAI_score', 'mean_ENOVA', 'VEOG', 'veog_corr', '-v7.3');
    fprintf('  Done.\n\n');

end

fprintf('===== Batch complete. %d file(s) processed. =====\n', nFiles);
