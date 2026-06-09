function dataset = preprocess_eeg(dataset)
%PREPROCESS_EEG Apply baseline correction and average re-referencing.

baselineMask = dataset.time < 0;
if ~any(baselineMask)
    error('Time vector must include a pre-stimulus baseline (< 0 s).');
end

baseline = mean(dataset.data(:, baselineMask, :), 2);
dataset.data = bsxfun(@minus, dataset.data, baseline);

channelMean = mean(dataset.data, 1);
dataset.data = bsxfun(@minus, dataset.data, channelMean);
end
