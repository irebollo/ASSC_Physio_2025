%% Tutorial Script 1: GLM_HERs
%
% This script loops over all subjects of the list 'included_subjects' and 
% computes a GLM to test to what extent a subjects' pre-stimulus heartbeat-evoked 
% response (HER) influenced their subsequent reaction times (RT) - more 
% specifically, here we use their residual RT (resRT) after regressing out 
% variance due to cardiac phase, inter-beat interval, and block number. Since 
% these are GLM residuals, they are already normally distributed and do not
% need to be normalized previous to our GLM. Of note, all data structures already 
% only contain clean valid trials of the conditions of interest.
%
% !!! 
% In practice, we only have the data for subject 26 (index 18 in our included_subjects 
% list) so we are commenting out the looping and just loading data and running
% the script for this subject.
% !!!
%
% INPUTS:
% The script loads, for each subject:
% - <filePrefix>_HERs.mat, containing data_HERs: a FieldTrip data structure
%       containing each trial's EEG epoch of interest (the last HER before 
%       stimulus onset)
% - <filePrefix>_resRTs.mat, containing data_resRTs: a 1 x nTrials vector
%       containing resRT for each trial
%
% FUNCTIONS CALLED:
% - TUTORIAL_setPaths: Go change the contents of this function to set paths
%       to the tutorial root directory and to your FieldTrip
% - FieldTrip
% - My_GLM
% 
% OUTPUTS:
% - <filePrefix>_GLM-results.mat, containing beta_HER_resRT: a FieldTrip
%       data structure containing, within the 'avg' field, this subjects'
%       beta timeseries resulting from the GLM. This beta_resRT expresses,
%       for each channel and each timepoint, to what extent the HER
%       co-varied with resRT.
% - <filePrefix>_group_HER.mat, containing HERs_allSubs: a 1 x nSubjects
%       cell array where each cell contains one subjects' HER averaged
%       across all trials in FieldTrip format
% - Figure of HER and beta timecourses
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


%% 

for subject = 18 % 1:length(included_subjects) % we skip the loop and only run on our example subject
    subject_nr = convertStringsToChars(included_subjects(subject));
    
    
    %% 1 - Load data
    
    load([paths.main 'Data' filesep 'sub-' subject_nr filesep 'sub-' subject_nr '_HERs.mat'], 'data_HERs');
    load([paths.main 'Data' filesep 'sub-' subject_nr filesep 'sub-' subject_nr '_resRTs.mat'], 'data_resRTs');
    
    
    %%% 1.1 Just to look at it: Compute average HER
    
    cfg = [];
    cfg.channel = {'all'};
    HER_avg = ft_timelockanalysis(cfg, data_HERs);
    
    cfg = [];
    cfg.layout = 'biosemi64.lay';
    cfg.interactive = 'yes';
    cfg.showoutline = 'yes';
    cfg.viewmode = 'butterfly';
    ft_multiplotER(cfg, HER_avg);
    title(['sub-' subject_nr])
    
    
    %%% 1.2 Prepare data for GLM
    
    % filling EEG data structure
    n_trials = length(data_HERs.trial);
    n_channels = length(data_HERs.label);
    n_timepoints = length(data_HERs.time{1,1});
    data_HERs_forGLM = nan(n_trials,n_channels,n_timepoints); % we are going to permute dimensions 2 and 3 below to fit requirements of My_glm.m
    for trial = 1:n_trials
        data_HERs_forGLM(trial,:,:) = data_HERs.trial{1,trial}(:,1:n_timepoints);
    end
    data_HERs_forGLM = permute(data_HERs_forGLM,[1,3,2]); % now data_eeg_forGLM = Trial x Time x Channel
    
    % if you have non-normally distributed regressors, normalize them here
    % if you need to create interaction terms, create them here
    
    
    %% 2 - Run GLM
    
    % GLM: HER ~ resRT
    % We try to explain variance in HERs by variance in (residual) reaction times
    
    regressors = [data_resRTs];
    zflag = 1; % zscore HER data across repetitions, separately per channel and per timepoint
    [beta, intercept, residuals] = My_glm(data_HERs_forGLM, regressors, zflag);
    
    % Store in structure
    beta_HER_resRT.avg       = beta; 
    beta_HER_resRT.dimord    = 'time_chan'; 
    beta_HER_resRT.intercept = intercept; 
    beta_HER_resRT.residuals = residuals;
    beta_HER_resRT.time      = data_HERs.time{1,1};
    beta_HER_resRT           = copyfields(data_HERs,beta_HER_resRT,{'label'});
    
    % Plot individual GLM result
    subplot(2,1,1)
    plot(HER_avg.time, HER_avg.avg)
    xlabel('Time (sec)')
    ylabel('Amplitude in channel (microV)')
    xlim([-0.03 0.6])
    xline(0,'-','Label','Rpeak')
    title('Avg HER')
    subplot(2,1,2)
    plot(beta_HER_resRT.time, beta_HER_resRT.avg)
    xlabel('Time (sec)')
    ylabel('beta(resRT)')
    xlim([-0.03 0.6])
    xline(0,'-','Label','Rpeak')
    title('beta(resRT)')
    sgtitle(['sub-' subject_nr ': GLM: HER ~ resRT'])
    print([paths.main 'Data' filesep 'group_results' filesep 'sub-' subject_nr '_GLM-results_beta'],'-dpng') 
    close
    
    % save individual GLM result
    save([paths.main 'Data' filesep 'group_results' filesep 'sub-' subject_nr '_GLM-results.mat'], 'beta_HER_resRT');
    clear beta_HER_resRT

    % Pool across subjects
    HERs_allSubs{subject} = HER_avg; clear HER_avg
end

% save([paths.main 'Data' filesep 'group_results' filesep 'group_HER.mat'], 'HERs_allSubs');
