% run_mmn_pipeline.m
% Minimal EEG/MMN analysis pipeline entrypoint.
%
% Update the paths and condition labels below, then run this script.

inputFile = 'data/example_eeg_epochs.mat';
outputFile = 'results/mmn_results.mat';
standardCondition = 'standard';
deviantCondition = 'deviant';

dataset = load_eeg_dataset(inputFile);
preprocessed = preprocess_eeg(dataset);
results = extract_mmn_erp(preprocessed, standardCondition, deviantCondition);
results.stats = cluster_permutation_test(results.mmn, results.time, 1000);

outputDir = fileparts(outputFile);
if ~isempty(outputDir) && ~isfolder(outputDir)
    mkdir(outputDir);
end

save(outputFile, 'results');
fprintf('Pipeline complete. Saved results to %s\n', outputFile);
