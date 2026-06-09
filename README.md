# MMN Gamified Phonetic Training
EEG/MMN analysis pipeline for a gamified phonetic learning study. Includes preprocessing, ERP extraction, and cluster-based permutation testing (FieldTrip) to examine neural discrimination of English vowel contrasts (/iː/–/ɪ/ and /æ/–/ɛ/) before and after training.
# Dependencies
MATLAB R2024b  
EEGLAB v2026  
FieldTrip fieldtrip-20250106  
GEDAI toolbox  (https://github.com/NeuroEngUAB/GEDAI)  
NoiseTools     (http://audition.ens.fr/adc/NoiseTools/)

# Data Acquisition
EEG system: NeuroScan (amplifier: SynAmps RT)  
Acquisition software: CURRY 8  
Sampling rate: 1000 Hz  
Channels: 30 EEG + 2 EOG (32-channel QuickCap, NeuroScan)  
Cap layout: standard 10-20 system  
# Trigger Codes and Channel Selection
Triggers were embedded in the EEG recording by the stimulus presentation software. The raw files contain 30 EEG channels, 2 EOG channels, and 1 trigger channel. The trigger channel is discarded during preprocessing and only the EEG and EOG channels are retained for analysis.

| Code   | Meaning               |
| ------ | --------------------- |
| 1      | Standard stimulus     |
| 2      | Deviant stimulus      |
| 800000 | Block boundary marker |

The 800000 marker is inserted by CURRY 8 at the start or boundary of each experimental block and is used by the pipeline to split the continuous EEG into separate blocks before epoching.

# Pipeline
01_filter_and_prepare_channels.m  
02_run_GEDAI_artefact_removal.m  
03_extract_ERPs_and_MMN.m  
04_MMN_by_contrast_pre_post.m  
05_cluster_permutation_MMN.m  
06_cluster_permutation_pre_vs_post.m

# Data Availability
Raw EEG data are not shared due to GDPR. Available on reasonable request.

# Acknowledgements
Analysis scripts were written in MATLAB and developed with the assistance of GitHub Copilot (Claude Sonnet 4.6), which was used for code generation, debugging, and pipeline development.

# Licence
MIT
