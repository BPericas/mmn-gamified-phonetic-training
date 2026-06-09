function results = extract_mmn_erp(dataset, standardLabel, deviantLabel)
%EXTRACT_MMN_ERP Compute ERP waveforms and MMN (deviant - standard).

conditionLabels = dataset.conditions;
standardIdx = conditionLabels == string(standardLabel);
deviantIdx = conditionLabels == string(deviantLabel);

if ~any(standardIdx)
    error('No trials found for standard condition "%s".', standardLabel);
end

if ~any(deviantIdx)
    error('No trials found for deviant condition "%s".', deviantLabel);
end

erpStandard = mean(dataset.data(:, :, standardIdx), 3);
erpDeviant = mean(dataset.data(:, :, deviantIdx), 3);

results = struct();
results.time = dataset.time;
results.erpStandard = erpStandard;
results.erpDeviant = erpDeviant;
results.mmn = erpDeviant - erpStandard;
results.mmnChannelAverage = mean(results.mmn, 1);
end
