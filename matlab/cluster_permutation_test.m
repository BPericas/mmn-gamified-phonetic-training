function stats = cluster_permutation_test(mmnByChannel, time, nPermutations, threshold)
%CLUSTER_PERMUTATION_TEST One-sample sign-flip cluster permutation test.
%
% Inputs:
%   mmnByChannel [channels x time]
%   time         [1 x time]
%   nPermutations number of permutations (default 1000)
%   threshold     cluster-forming |t| threshold (default 2.0)
%
% Note: channels are treated as exchangeable observations for this
% minimal implementation.

if nargin < 3
    nPermutations = 1000;
end
if nargin < 4
    threshold = 2.0;
end

nChannels = size(mmnByChannel, 1);

sampleMean = mean(mmnByChannel, 1);
sampleStd = std(mmnByChannel, 0, 1);
tValues = sampleMean ./ (sampleStd ./ sqrt(nChannels) + eps);

clusters = find_clusters(abs(tValues) > threshold);
clusterMasses = zeros(1, numel(clusters));
for i = 1:numel(clusters)
    clusterMasses(i) = sum(abs(tValues(clusters{i})));
end

nullDistribution = zeros(1, nPermutations);
for p = 1:nPermutations
    flips = sign(randn(nChannels, 1));
    permuted = mmnByChannel .* flips;
    pMean = mean(permuted, 1);
    pStd = std(permuted, 0, 1);
    pT = pMean ./ (pStd ./ sqrt(nChannels) + eps);

    pClusters = find_clusters(abs(pT) > threshold);
    if isempty(pClusters)
        nullDistribution(p) = 0;
    else
        pMasses = cellfun(@(idx) sum(abs(pT(idx))), pClusters);
        nullDistribution(p) = max(pMasses);
    end
end

pValues = ones(1, numel(clusterMasses));
for i = 1:numel(clusterMasses)
    pValues(i) = mean(nullDistribution >= clusterMasses(i));
end

stats = struct();
stats.time = time;
stats.tValues = tValues;
stats.threshold = threshold;
stats.clusters = clusters;
stats.clusterMasses = clusterMasses;
stats.clusterPValues = pValues;
stats.nullDistribution = nullDistribution;
end

function clusters = find_clusters(mask)
if ~any(mask)
    clusters = {};
    return;
end

idx = find(mask);
edges = [1, find(diff(idx) > 1) + 1, numel(idx) + 1];
nClusters = numel(edges) - 1;
clusters = cell(1, nClusters);
for clusterIdx = 1:numel(edges) - 1
    range = edges(clusterIdx):edges(clusterIdx + 1) - 1;
    clusters{clusterIdx} = idx(range);
end
end
