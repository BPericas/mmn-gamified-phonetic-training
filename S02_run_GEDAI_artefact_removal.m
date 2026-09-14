%% S02_run_GEDAI_artefact_removal.m
%
%  Second preprocessing step in the MMN analysis pipeline.
%
%  Runs GEDAI (Generalised Eigenvalue Decomposition for Artefact
%  Identification) on all filtered files produced by script 01. ICA is then
%  run on the GEDAI-cleaned data and ICLabel is used to identify and remove
%  eye-blink components, after which a robust average rereference is
%  applied using NoiseTools (nt_rereference).
%
%  Input:   <subjectID>_filtnotrig.set   (output of S01_filter_and_prepare_channels.m)
%           <subjectID>_VEOG.mat         (VEOG time course saved by script 01)
%  Output:  <subjectID>_GEDAIclean.set   — artefact-cleaned, ICA blink-removed,
%                                          rereferenced EEG
%           <subjectID>_GEDAIclean.mat   — MATLAB snapshot including GEDAI
%                                          quality metrics, ICLabel classifications,
%                                          removed blink components and VEOG correlations
%
%  Dependencies
%  ------------
%  - GEDAI toolbox        (https://github.com/NeuroEngUAB/GEDAI)
%  - NoiseTools           (http://audition.ens.fr/adc/NoiseTools/)
%  - EEGLAB ICLabel plugin (eye-component classification)
%
%  Notes
%  -----
%  - Eye-blink components are identified using ICLabel (a neural network
%    classifier built into EEGLAB) run on the GEDAI-cleaned data — no VEOG
%    channel is required for removal. Components are only removed if they
%    are classified as "Eye" with probability >= eye_prob_threshold AND
%    have a frontopolar (Fp1/Fp2) topography loading >= fp_threshold times
%    the across-channel mean, which filters out ICLabel false positives.
%  - VEOG is not used for artefact removal. It is loaded solely to compute
%    a post-ICA residual correlation check, saved for inspection — no
%    automatic exclusion is applied based on it.
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
%% ICA / ICLabel eye-removal parameters — adjust if needed

eye_prob_threshold = 0.80;   % remove ICs classified as eye with >= 80% probability
fp_labels          = {'Fp1','Fp2','FP1','FP2'};   % frontopolar channel labels to check
fp_threshold       = 2.0;    % Fp loading must be >= 2x the across-channel mean

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
    % Diagnostic only — used for the post-ICA residual correlation check,
    % not for artefact removal (ICLabel handles that without external channels).
    veogFile = fullfile(fileDir, [subjectID '_VEOG.mat']);
    if exist(veogFile, 'file')
        tmp     = load(veogFile, 'VEOG');
        VEOG    = tmp.VEOG(:);
        hasVEOG = true;
        fprintf('  VEOG loaded (%d samples) — will be used for post-ICA QC only.\n', length(VEOG));
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

    %% STEP 4: ICA blink removal using ICLabel
    %
    %  Runs ICA on the GEDAI-cleaned data, then uses ICLabel to identify eye
    %  components by their topography and activation characteristics. ICLabel
    %  assigns probability scores across 7 classes; components with eye
    %  probability >= eye_prob_threshold AND a frontopolar topography are removed.
    %
    %  ICLabel class order: Brain | Muscle | Eye | Heart | LineNoise | ChanNoise | Other
    blink_ics = [];   % IC indices removed
    ic_labels = [];   % full ICLabel probability matrix [nICs x 7]

    fprintf('  Running ICA (runica, extended)...\n');
    EEG = pop_runica(EEG, 'icatype', 'runica', 'extended', 1);
    EEG = eeg_checkset(EEG);

    fprintf('  Running ICLabel classifier...\n');
    EEG = pop_iclabel(EEG, 'default');

    % ic_labels rows = ICs, columns = [Brain Muscle Eye Heart LineNoise ChanNoise Other]
    ic_labels = EEG.etc.ic_classification.ICLabel.classifications;

    % Flag eye components by ICLabel probability (column 3 = Eye probability)
    EEG = pop_icflag(EEG, ...
        [NaN NaN;            ...   % Brain      — never reject
         NaN NaN;            ...   % Muscle     — not targeting here
         eye_prob_threshold 1; ... % Eye        — flag if >= threshold
         NaN NaN;            ...   % Heart
         NaN NaN;            ...   % Line Noise
         NaN NaN;            ...   % Channel Noise
         NaN NaN]);                % Other

    eye_ics_iclabel = find(EEG.reject.gcompreject);

    % Additional filter: require frontopolar topography, since blink ICs
    % project maximally to Fp1/Fp2. This removes ICLabel false positives
    % (e.g. muscle or brain ICs with coincidentally high eye probability).
    fp_chans = find(ismember({EEG.chanlocs.labels}, fp_labels));
    ic_topo  = abs(EEG.icawinv);   % [channels x ICs]

    if ~isempty(fp_chans) && ~isempty(eye_ics_iclabel)
        fp_weight = mean(ic_topo(fp_chans, :), 1)' ./ (mean(ic_topo, 1)' + eps);
        blink_ics = eye_ics_iclabel(fp_weight(eye_ics_iclabel) >= fp_threshold);
    elseif isempty(fp_chans)
        warning('Fp1/Fp2 not found in channel labels — using ICLabel probability only.');
        blink_ics = eye_ics_iclabel;
    else
        blink_ics = [];
    end

    fprintf('  ICLabel eye ICs (prob >= %.2f)        : %s\n', ...
        eye_prob_threshold, num2str(eye_ics_iclabel'));
    fprintf('  Blink ICs (ICLabel + frontopolar topo): %s\n', num2str(blink_ics'));

    if ~isempty(blink_ics)
        EEG = pop_subcomp(EEG, blink_ics, 0);
        EEG = eeg_checkset(EEG);
        fprintf('  Removed %d eye IC(s).\n', numel(blink_ics));
    else
        fprintf('  No eye ICs above threshold — data unchanged.\n');
    end

    %% STEP 5: Post-ICA VEOG residual correlation check
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
        fprintf('  Max |VEOG corr| post-ICA: %.3f (channel %s)\n', ...
            max_r, EEG.chanlocs(max_ch).labels);
        if max_r > 0.20
            warning('  Residual VEOG correlation high (%.3f) for %s — inspect IC topographies manually.', ...
                max_r, subjectID);
        end
    else
        veog_corr = [];
    end

    %% STEP 6: Robust average rereference (NoiseTools)
    %
    %  nt_rereference computes a weighted average reference that is robust
    %  to channels with high noise, unlike a simple average reference.
    fprintf('  Applying robust average reference (NoiseTools)...\n');
    x        = double(EEG.data)';   % NoiseTools expects [samples x channels]
    x        = nt_rereference(x);
    EEG.data = x';                  % back to [channels x samples]
    EEG      = eeg_checkset(EEG);

    %% STEP 7: Save cleaned data
    fprintf('  Saving %s_GEDAIclean ...\n', subjectID);
    EEG = pop_saveset(EEG, ...
        'filename', [subjectID '_GEDAIclean.set'], ...
        'filepath', fileDir);

    save(outMatFile, 'EEG', 'subjectID', 'SENSAI_score', 'mean_ENOVA', 'VEOG', ...
         'veog_corr', 'blink_ics', 'ic_labels', '-v7.3');
    fprintf('  Done.\n\n');

end

fprintf('===== Batch complete. %d file(s) processed. =====\n', nFiles);

