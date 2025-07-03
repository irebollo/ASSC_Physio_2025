%% Tutorial Script 2: GroupStats_HERs
%
% This script performs group statistics to test whether pre-stimulus heartbeat-evoked 
% potentials (HERs) significantly predicted subsequent reaction times (resRT). 
% Specifically, it performs permutation-based cluster statistics on subjects' beta
% space-time-series resulting from individual GLM's (output of Analyse_HERs) 
% using ft_clusterstats. It loops over subjects to load everyones' 
% GLM results, pools them into a common structure, and then performs group stats. 
% Group results are plotted. 
%
% INPUTS:
% The script loads, for each subject:
% - <filePrefix>_subjectNumber_GLM-results.mat, containing beta_HER_resRT: a FieldTrip
%       data structure containing, within the 'avg' field, this subjects'
%       beta timeseries resulting from the GLM run in Analyse_HERs. This 
%       beta_resRT expresses, for each channel and each timepoint, to what 
%       extent the HER co-varied with resRT.
%
% FUNCTIONS CALLED:
% - TUTORIAL_setPaths: Go change the contents of this function to set paths
%       to the tutorial root directory and to your FieldTrip
% - FieldTrip
% - get_clusterInfo
%
% OUTPUTS:
% - <filePrefix>_GLM-results.mat, containing stat_resRT: Results of group 
%       statistics run on beta timeseries of all subjects. Output of ft_clusterstats.
% - Plots of topographies and beta timecourses for all significant
%       clusters, indicating where and when HER significantly co-varied with
%       resRT
%
% AUTHOR: Marie Loescher - July 2025


%% INIT

clear all

% Defining and adding paths
paths = TUTORIAL_setPaths();

% Fieldtrip
restoredefaultpath;
addpath([paths.fieldtrip])
ft_defaults

addpath(genpath([paths.main]));

load([paths.main 'Data' filesep 'included_subjects.mat'],'included_subjects')
outputPath = [paths.main 'Data' filesep 'group_results' filesep 'group_GLM-results'];


%% 1 - Pool data of all subjects

for subject = 1:length(included_subjects) 
    subject_nr = convertStringsToChars(included_subjects(subject));

    load([paths.main 'Data' filesep 'group_results' filesep 'sub-' subject_nr '_GLM-results.mat'], 'beta_HER_resRT')

    beta_HER_resRT_all{subject}.avg     = transpose(beta_HER_resRT.avg); 
    beta_HER_resRT_all{subject}.time    = beta_HER_resRT.time; 
    beta_HER_resRT_all{subject}.label   = beta_HER_resRT.label; 
    beta_HER_resRT_all{subject}.dimord  = 'chan_time'; 

    % creating matching data that contains zero to perform the statistical analysis against
    beta_ZEROES{subject}.avg     = transpose(zeros(size(beta_HER_resRT.avg))); 
    beta_ZEROES{subject}.time    = beta_HER_resRT.time; 
    beta_ZEROES{subject}.label   = beta_HER_resRT.label; 
    beta_ZEROES{subject}.dimord  = 'chan_time'; 

    clear beta_HER_resRT
end


%% 2 - Perform group statistics: Compare betas to 0

%%% prepare neigbhours
cfg             = [];
cfg.channel     = 'all';
cfg.method      = 'template';
cfg.template    = 'biosemi64_neighb.mat';
cfg.layout      = 'biosemi64.lay';  
neighbours      = ft_prepare_neighbours(cfg);

% create design matrix
n_subjects = size(included_subjects,2);
design = zeros(2,2*n_subjects);
for i = 1:n_subjects
    design(1,i) = i;
end
for i = 1:n_subjects
    design(1,n_subjects+i) = i;
end
design(2,1:n_subjects) = 1;
design(2,n_subjects+1:2*n_subjects) = 2;

%%% Testing against zeroes
cfg = [];
cfg.channel = {'all'};
cfg.latency = [0.2 0.6]; % R peak = 0
cfg.method = 'montecarlo';
cfg.statistic = 'ft_statfun_depsamplesT';
cfg.avgovertime = 'no';
cfg.correctm = 'cluster';
cfg.parameter = 'avg';
cfg.clusterstatistic = 'maxsum';
cfg.tail = 0; 
cfg.clustertail = 0; 
cfg.alpha = 0.05;
cfg.correcttail = 'alpha'; % because we are performing two-sided t-tests
cfg.numrandomization = 2000;
cfg.design = design;
cfg.neighbours = neighbours;
cfg.uvar  = 1; % the 1st row in cfg.design contains the subject number
cfg.ivar  = 2; % the 2nd row in cfg.design contains the independent variable
cfg.spmversion = 'spm12'; % MEX files in SPM 8 not compatible, must use spm12
cfg.clusteralpha = 0.01;
cfg.minnbchan = 1;

stat_resRT  = ft_timelockstatistics(cfg, beta_HER_resRT_all{:}, beta_ZEROES{:});

save([outputPath '.mat'], 'stat_resRT');


%% 3 - Plotting group results

load([paths.main 'Data' filesep 'group_results' filesep 'group_HER.mat'], 'HERs_allSubs'); % load all subjects' HERs for plot

% Retrieve info on significant timepoints and channels of the positive clusters for plotting
[timeStart_resRT_1, timeEnd_resRT_1, channels_i_resRT_1, maxChannel_i_resRT_1] = get_clusterInfo(stat_resRT, 'pos', 1); % first cluster
[timeStart_resRT_2, timeEnd_resRT_2, channels_i_resRT_2, maxChannel_i_resRT_2] = get_clusterInfo(stat_resRT, 'pos', 2); % second cluster
channels_i_resRT_bothClusters = intersect(channels_i_resRT_1,channels_i_resRT_2); % these are the channels that are significant in both clusters


%%% Topoplot of cluster 1

cfg                         = []; 
cfg.parameter               = 'stat';
cfg.layout                  = 'biosemi64.lay';
cfg.colorbar                = 'yes';
cfg.colormap                = ft_colormap('RdBu');
cfg.zlim                    = [-4 4];
cfg.highlight               = 'on';
cfg.highlightsymbol         = 'o';
cfg.highlightcolor          = 'white';
cfg.highlightchannel        = channels_i_resRT_1;
cfg.highlightsize           = 9;
cfg.markercolor             = 'none';
cfg.visible                 = 'on';
cfg.xlim                     = [timeStart_resRT_1 timeEnd_resRT_1];
cfg.colorbar = 'SouthOutside';
ft_topoplotER(cfg,stat_resRT)
set(gcf, 'Color', 'w')
print([outputPath '_topo_1'],'-dpng') 
close


%%% Topoplot of cluster 2

cfg                         = []; 
cfg.parameter               = 'stat';
cfg.layout                  = 'biosemi64.lay';
cfg.colorbar                = 'yes';
cfg.colormap                = ft_colormap('RdBu');
cfg.zlim                    = [-4 4];
cfg.highlight               = 'on';
cfg.highlightsymbol         = 'o';
cfg.highlightcolor          = 'white';
cfg.highlightchannel        = channels_i_resRT_2;
cfg.highlightsize           = 9;
cfg.markercolor             = 'none';
cfg.visible                 = 'on';
cfg.xlim                     = [timeStart_resRT_2 timeEnd_resRT_2];
cfg.colorbar = 'SouthOutside';
ft_topoplotER(cfg,stat_resRT)
set(gcf, 'Color', 'w')
print([outputPath '_topo_2'],'-dpng') 
close


%%% Beta timecourse over significant electrodes common in both clusters

subplot(2,1,1)

% Compute mean beta over channels per subject
for subject = 1:size(beta_HER_resRT_all,2) 
    % do this to copy all fields first
    chanAvg_1{subject} = beta_HER_resRT_all{subject};
    % and now replace the avg field with data from only channels of interest
    chanAvg_1{subject}.avg = mean(beta_HER_resRT_all{subject}.avg(channels_i_resRT_bothClusters,:));
end
% Then compute Grand average (GA) and Standard error (SE) of betas over subjects
cfg = [];
cfg.channel = 'all'; % this data now already contains only average over channels of interest
cfg.latency = 'all';
cfg.parameter = 'avg';
GA_chanAvg_1 = ft_timelockgrandaverage(cfg, chanAvg_1{:});
SE_chanAvg_1 = (sqrt(GA_chanAvg_1.var(:))/sqrt((GA_chanAvg_1.dof(1))))';

% Compute bounds of SE plot
time = GA_chanAvg_1.time';
low_1 = (GA_chanAvg_1.avg - SE_chanAvg_1)'; % lower bound of SE for plot
high_1 = (GA_chanAvg_1.avg + SE_chanAvg_1)'; % upper bound of SE for plot

% plot standard error area
SE_plot_1 = patch([time; time(end:-1:1); time(1)], [low_1; high_1(end:-1:1); low_1(1)], 'k');
set(SE_plot_1,'facecolor','#34689F','edgecolor','none','FaceAlpha',0.5);
hold on
% plot channel average
line(time, GA_chanAvg_1.avg, 'Color', '#34689F', 'linewidth', 1);
% plot significant timepoints
line([timeStart_resRT_1,timeEnd_resRT_1],[-0.01, -0.01], 'Color','black', 'LineWidth',2)
line([timeStart_resRT_2,timeEnd_resRT_2],[-0.01, -0.01], 'Color','black', 'LineWidth',2)
% Add transparency to un-tested timewindow
lowlight = patch([-0.03 0.2 0.2 -0.03], [-1 -1 1 1], 'white');
set(lowlight, 'edgecolor','none','FaceAlpha',0.5);
% figure settings
ylim([-0.08 0.08])
xlim([-0.03 0.6]);
ylabel('beta(resRT)')
xlabel('Time (sec)')
xline(0,'-','linewidth',0.5,'Label','Rpeak')
yline(0,'--','linewidth',0.5)
title('Group GLM result: beta(resRT) of GLM: HER ~ resRT over significant channels')


%%% HER timecourse over significant electrodes common in both clusters

subplot(2,1,2)

% Compute mean HER over channels per subject
for subject = 1:size(HERs_allSubs,2) 
    % do this to copy all fields first
    chanAvg_1{subject} = HERs_allSubs{subject};
    % and now replace the avg field with data from only channels of interest
    chanAvg_1{subject}.avg = mean(HERs_allSubs{subject}.avg(channels_i_resRT_bothClusters,:));
end
% Then compute Grand average (GA) and Standard error (SE) of HERs over subjects
cfg = [];
cfg.channel = 'all'; % this data now already contains only average over channels of interest
cfg.latency = 'all';
cfg.parameter = 'avg';
GA_chanAvg_1 = ft_timelockgrandaverage(cfg, chanAvg_1{:});
SE_chanAvg_1 = (sqrt(GA_chanAvg_1.var(:))/sqrt((GA_chanAvg_1.dof(1))))';

% Compute bounds of SE plot
time = GA_chanAvg_1.time';
low_1 = (GA_chanAvg_1.avg - SE_chanAvg_1)'; % lower bound of SE for plot
high_1 = (GA_chanAvg_1.avg + SE_chanAvg_1)'; % upper bound of SE for plot

% plot standard error area
SE_plot_1 = patch([time; time(end:-1:1); time(1)], [low_1; high_1(end:-1:1); low_1(1)], 'k');
set(SE_plot_1,'facecolor','#34689F','edgecolor','none','FaceAlpha',0.5);
hold on
% plot channel average
line(time, GA_chanAvg_1.avg, 'Color', '#34689F', 'linewidth', 1);
% plot significant timepoints
line([timeStart_resRT_1,timeEnd_resRT_1],[-0.01, -0.01], 'Color','black', 'LineWidth',2)
line([timeStart_resRT_2,timeEnd_resRT_2],[-0.01, -0.01], 'Color','black', 'LineWidth',2)
% Add transparency to un-tested timewindow
lowlight = patch([-0.03 0.2 0.2 -0.03], [-1 -1 1 1], 'white');
set(lowlight, 'edgecolor','none','FaceAlpha',0.5);
% figure settings
ylim([-0.6 0.6])
xlim([-0.03 0.6]);
ylabel('HER amplitude (microV)')
xlabel('Time (sec)')
xline(0,'-','linewidth',0.5,'Label','Rpeak')
yline(0,'--','linewidth',0.5)
title('Average HER over significant channels')

%%% save
print([outputPath '_beta'],'-dpng') 
close

