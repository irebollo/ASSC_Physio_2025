function [maxstat_pos, maxstat_neg] = INTESPACE_find_maxstat(clusterstat_result)

% This function finds largest cluster statistic (for instance maxSumT) for
% positive and negative candidate clusters contained in clusterstat_result, output of
% ft_clusterstats. If there is no candidate cluster, cluster statistic is
% set to 0.

% INPUTS:
% - clusterstat_result: output of ft_clusterstats

% OUTPUTS:
% - maxstat_pos and maxstat_neg, value of largest cluster statistic found
%       for positive and negative candidate clusters in clusterstat_result,
%       respectively

% AUTHOR:
% Marie Loescher - May 2024


%%

if isfield(clusterstat_result,'posclusters')
    if size(clusterstat_result.posclusters,2) ~= 0
        maxstat_pos = clusterstat_result.posclusters(1).clusterstat; % the first is always the largest
    else
        maxstat_pos = 0; % attribute 0 value if no positive cluster is found
    end
else
    maxstat_pos = 0;
end

if isfield(clusterstat_result,'negclusters')
    if size(clusterstat_result.negclusters,2) ~= 0
        maxstat_neg = clusterstat_result.negclusters(1).clusterstat;
    else
        maxstat_neg = 0;
    end
else
    maxstat_neg = 0;
end




end