function dataset = load_eeg_dataset(filePath)
%LOAD_EEG_DATASET Load EEG epoch data from a MAT file.
% Expected variables in file:
%   data       [channels x time x trials]
%   time       [1 x time] in seconds
%   fs         sampling frequency in Hz
%   conditions [1 x trials] cell array/string array of condition labels

loaded = load(filePath);
required = {'data', 'time', 'fs', 'conditions'};

for i = 1:numel(required)
    fieldName = required{i};
    if ~isfield(loaded, fieldName)
        error('Missing required variable "%s" in %s', fieldName, filePath);
    end
end

dataset = struct();
dataset.data = loaded.data;
dataset.time = loaded.time(:)';
dataset.fs = loaded.fs;
dataset.conditions = string(loaded.conditions);
end
