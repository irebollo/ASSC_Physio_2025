%% Tutorial Preparation Script
% To select data from INTESPACE

clear all

% Defining and adding paths
paths = INTESPACE_setPaths();

restoredefaultpath;
addpath([paths.main 'Scripts' filesep 'Toolboxes' filesep 'fieldtrip-20221212' filesep])
ft_defaults

addpath(genpath([paths.main 'Data' filesep]));
addpath(genpath([paths.main 'Scripts' filesep 'Preprocessing']));
addpath(genpath([paths.main 'Scripts' filesep 'Analysis']));
addpath(genpath([paths.main 'Scripts' filesep 'Toolboxes' filesep 'Other']));
addpath(genpath([paths.main 'Scripts' filesep 'Toolboxes' filesep 'brainstorm3-master']));

path2tutorial = '\\129.199.81.50\intespace\4TUTORIAL\';
load([paths.main 'Data' filesep 'included_subjects.mat'],'included_subjects')

% Set up the blocks ran for each subject
for subject = 1:length(included_subjects) 
    subject_nr = convertStringsToChars(included_subjects(subject));
    if strcmp(subject_nr,"13") 
        subject_blocks{subject} = ["01","02","03","04", "04-2","05","RS"]; % block recording was interrupted into two parts
    elseif strcmp(subject_nr,"47") 
        subject_blocks{subject} = ["01","02","02-2","03", "04","05","RS"]; % block recording was interrupted into two parts
    else
        subject_blocks{subject} = ["01","02","03","04","05","RS"];
    end
end


for subject = 1:length(included_subjects) 
    subject_nr = convertStringsToChars(included_subjects(subject))
    [~, subject_paths] = INTESPACE_setPaths(subject_nr);
    all_blocks = subject_blocks{subject};
    

    %% Get data from INTESPACE

%     load([paths.preproc 'sub-' subject_nr filesep 'Behaviour' filesep 'sub-' subject_nr '_desc-preproc_beh.mat'], 'data_beh_OK');
%     distanceConditions = data_beh_OK(:,4); clear data_beh_OK 
%     load([paths.results 'beh_ECG_LMMs_GLMs' filesep 'sub-' subject_nr '_GLM-RTvTAT_beta.mat'], 'resRT')
%     load([subject_paths.preproc.eeg 'sub-' subject_nr '_HER-R600_bp-0525_cond-TAT_eeg.mat'], 'data_eeg_bp_Rpeaks_TATclean')
%     load([paths.results 'EEG_HER_GLMs' filesep 'sub-' subject_nr '_GLM-HERv4-Dist_R-last600_bp-0525_cond-TAT_beta.mat'], 'trials_of_HERs_UsedInGLM')
%     load([paths.preproc 'sub-' subject_nr filesep 'sub-' subject_nr '_shuffledHER-600_cond-TAT_events.mat'], 'trialDef_Rpeaks_allShuffles_allBlocks');
%     load([paths.preproc 'sub-' subject_nr filesep 'sub-' subject_nr '_HER-600_cond-TAT_events.mat'], 'trialDef_Rpeaks_allBlocks');

    for b = 1:length(all_blocks)-1 % skip RS block
        block_nr = convertStringsToChars(all_blocks(b));
        load([subject_paths.preproc.eeg filesep 'sub-' subject_nr '_block-' block_nr '_desc-ICAcardio_bp-0525_eeg.mat'],'data_eeg_ICAcardio_bp');
        save([path2tutorial 'Data' filesep 'sub-' subject_nr filesep 'sub-' subject_nr '_block-' block_nr '_desc-ICAcardio_bp-0525_eeg.mat'],'data_eeg_ICAcardio_bp');
    end


%     %% Organize data for tutorial
% 
%     % select only trialDefs of interest
%     for b = 1:length(all_blocks)-1 % skip RS block
%         block_nr = convertStringsToChars(all_blocks(b));
%         if regexp(['sub-' subject_nr '_block-' block_nr], regexptranslate('wildcard', '*sub-16_block-03*')) % this block was faulty. Leads to error because no valid trials left so skip it here already
%             continue
%         else
%             trialDef_Rpeaks{1,b} = trialDef_Rpeaks_allBlocks{1,b}((trialDef_Rpeaks_allBlocks{1,b}(:,5) == -1) & (trialDef_Rpeaks_allBlocks{1,b}(:,7) > 0),:);
%         end
%     end
%     trialDef_Rpeaks_allBlocks = trialDef_Rpeaks; clear trialDef_Rpeaks
% 
%     % select only HERs of interest
%     cfg = [];
%     cfg.trials = (ismember(data_eeg_bp_Rpeaks_TATclean.trialinfo(:,1),trials_of_HERs_UsedInGLM)) & (data_eeg_bp_Rpeaks_TATclean.trialinfo(:,2)==-1);
%     data_HERs = ft_selectdata(cfg,data_eeg_bp_Rpeaks_TATclean);
%     
%     % select only resRT and condition of trials of interest
%     for e = 1:length(trials_of_HERs_UsedInGLM)
%         resRT_perHER(e,1) = resRT(trials_of_HERs_UsedInGLM(e)); % --> one RT value per HER epoch
%         cond_perHER(e,1) = distanceConditions(trials_of_HERs_UsedInGLM(e)); % --> one distance code per HER epoch
%     end
% 
%     % rename for clarity
%     data_resRTs = resRT_perHER;
%     trialinfo.trialindex = trials_of_HERs_UsedInGLM;
%     trialinfo.soundDistance = cond_perHER;
%     trialDef_Rpeaks_allPermutations_allBlocks = trialDef_Rpeaks_allShuffles_allBlocks;
% 
%     % Save in tutorial folder
%     if ~exist([path2data 'Data' filesep 'sub-' subject_nr filesep])
%         mkdir([path2data 'Data' filesep 'sub-' subject_nr filesep]);
%     end
%     save([path2data 'Data' filesep 'sub-' subject_nr filesep 'sub-' subject_nr '_HERs.mat'], 'data_HERs');
%     save([path2data 'Data' filesep 'sub-' subject_nr filesep 'sub-' subject_nr '_trialinfo.mat'], 'trialinfo');
%     save([path2data 'Data' filesep 'sub-' subject_nr filesep 'sub-' subject_nr '_resRTs.mat'], 'data_resRTs');
%     save([path2data 'Data' filesep 'sub-' subject_nr filesep 'sub-' subject_nr '_HER_event.mat'], 'trialDef_Rpeaks_allBlocks');
%     save([path2data 'Data' filesep 'sub-' subject_nr filesep 'sub-' subject_nr '_permutedHER_events.mat.mat'], 'trialDef_Rpeaks_allPermutations_allBlocks');
% 
%     % Compute HER
%     cfg = [];
%     cfg.channel = {'all'};
%     HER_avg = ft_timelockanalysis(cfg, data_HERs);
%     HERs_allSubs{subject} = HER_avg;
% 
%     clear e subject_paths trials_of_HERs_UsedInGLM data_eeg_bp_Rpeaks_TATclean resRT distanceConditions resRT_perHER cond_perHER HER_avg trialDef_Rpeaks_allBlocks
end

save([path2tutorial 'Data' filesep 'group_results' filesep 'group_HER.mat'], 'HERs_allSubs');
