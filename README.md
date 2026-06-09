# mmn-gamified-phonetic-training
EEG/MMN analysis pipeline for a gamified phonetic learning study. Includes preprocessing, ERP extraction, and cluster-based permutation testing to examine neural discrimination of English vowel contrasts (/iː/–/ɪ/ and /æ/–/ɛ/) before and after training.

## MATLAB pipeline

This repository now includes a freely available MATLAB EEG/MMN workflow under `matlab/`:

- `run_mmn_pipeline.m` – entrypoint script for the full pipeline
- `load_eeg_dataset.m` – loads epoch data from a `.mat` file
- `preprocess_eeg.m` – baseline correction and average re-reference
- `extract_mmn_erp.m` – computes standard/deviant ERPs and MMN
- `cluster_permutation_test.m` – cluster-based sign-flip permutation statistics

### Expected input format

Provide a `.mat` file with these variables:

- `data`: EEG epochs with shape `[channels x time x trials]`
- `time`: time vector in seconds (must include a pre-stimulus interval `< 0`)
- `fs`: sampling frequency (Hz)
- `conditions`: condition label per trial (for example `standard` / `deviant`)

### Running

1. Edit paths/labels in `run_mmn_pipeline.m`.
2. Run the script in MATLAB.
3. Results are saved to `results/mmn_results.mat`.
