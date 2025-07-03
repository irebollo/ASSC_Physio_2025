function [timeStart, timeEnd, channels_i, maxChannel_i] = get_clusterInfo(stat, clusterType, clusterNumber)

% This function gives the start and end times and channels of the
% significant cluster specified by clusterType and clusterNumber.

% INPUTS:
% - stat: output of ft_clusterstats
% - clusterType: 'pos' or 'neg'
% - clusterNumber: A number; the clusters are ordered from largest to smallest in stat

% OUTPUT:
% - timeStart and timeEnd: time of start and end of cluster, as retrieved from the time field of stat
% - channels_i: all channels that are part of the significant cluster
% - maxChannel_i: channel with the highest cumulated test statistic over
%       significant timepoints

% AUTHOR: 
% Marie Loescher - Sept 2024

%%

if strcmp(clusterType, 'pos')
    thisClusterSamples = stat.posclusterslabelmat == clusterNumber;
    p = stat.posclusters(clusterNumber).prob;
elseif strcmp(clusterType, 'neg')
    thisClusterSamples = stat.negclusterslabelmat == clusterNumber;
    p = stat.negclusters(clusterNumber).prob;
else
    error('You should provide ''pos'' or ''neg'' as an argument for clusterType')
end

[~,sampleStart] = find(thisClusterSamples == 1,1,'first');
[~,sampleEnd] = find(thisClusterSamples == 1,1,'last');
timeStart = stat.time(sampleStart);
timeEnd = stat.time(sampleEnd);
channels_i = find(sum(thisClusterSamples,2)>0) ; % these are all channels with any significant timepoints
for c = 1:length(channels_i)
    c_i = channels_i(c);
    maxSumT(c) = sum(stat.stat(c_i,sampleStart:sampleEnd));
end
if strcmp(clusterType, 'pos')
    [~,max_i] = max(maxSumT);
elseif strcmp(clusterType, 'neg')
    [~,max_i] = min(maxSumT);
end
maxChannel_i = channels_i(max_i);

disp([clusterType ' cluster ' num2str(clusterNumber) ': p = ' num2str(p) ', t = ' num2str(timeStart) ' to ' num2str(timeEnd)])


end