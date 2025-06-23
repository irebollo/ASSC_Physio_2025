function [EGG_phaseXVolume] = compute_phaseXvolume(EGG_filtered, MRI_TR)

%%% Edited for Brain-Body Waves 2024 LB 

%% Downsample filtered EGG to MRI sampling rate 

disp('Resampling...')
cfg = [];  
cfg.detrend = 'no'; 
cfg.demean = 'yes';
cfg.resamplefs= 1/MRI_TR; % e.g., one volume every 1.4 seconds (TR) (1/1.4 = 0.7143 Hz)
EGG_filt_downsampled_MRI = ft_resampledata(cfg,EGG_filtered);

%% compute phase of the gastric cycle and the amplitude envelope of the filtered EGG, using the Hilbert method.

% Hilbert-transform of the EGG
EGG_phase = EGG_filt_downsampled_MRI; % just for EGG structure
% phase
EGG_phase.trial{1}      = angle(hilbert(EGG_filt_downsampled_MRI.trial{1}'))'; %plot(EGG_phase.time{1,:}, EGG_phase.trial{1,1}(1,:))
% amplitude envelope
EGG_phase.trial{1}(2,:) = abs(hilbert(EGG_filt_downsampled_MRI.trial{1}'))'; %plot(EGG_phase.time{1,:}, EGG_phase.trial{1,1}(2,:))   
% downsampled filtered EGG 
EGG_phase.trial{1}(4,:) = EGG_filt_downsampled_MRI.trial{1,1}(1,:);
% phase degrees
EGG_phase.trial{1}(3,:) = radtodeg(EGG_phase.trial{1}(1,:));
% labels
EGG_phase.label         = {'phase', 'amplitude', 'filteredEGG', 'phase_degrees'};


%% Get average phase value per volume

nVolumes = (length(EGG_filt_downsampled_MRI.trial{1,1})/EGG_filt_downsampled_MRI.fsample) / MRI_TR;
if round(nVolumes) == length(EGG_phase.trial{1})
    EGG_phaseXVolume= EGG_phase; 
end

% % If have TR markers/triggers in EGG recording:
% cfg = []; % configuration structure 
% cfg.dataset = subj_filename; % name of EGG file 
% cfg.trialfun = 'ft_trialfun_general'; % function
% cfg.trialdef.eventtype = 'Response'; % type of trigger
% cfg.trialdef.eventvalue = 'R128'; % name of TR trigger
% cfg.trialdef.prestim = 1; % prestim win seconds for epoching
% cfg.trialdef.poststim = 1; % poststim win seconds for epoching
% markersInFieldtrip = ft_definetrial(cfg); % load TR triggers
% 
% % downsample markers/triggers
% markers_downsampled_MRI = prepro_egg_downsampleMarkers(EGG_raw,EGG_downsampled_MRI,markersInFieldtrip); % Downsample markers to new MRI SR
% 
% nVolumes=length(markers_downsampled_MRI);
% 
% % compute average EGG phase per MRI volume
% EGG_phaseXVolume = zeros(1,nVolumes);
% for iTrial=1:nVolumes
%     EGG_phaseXVolume(1,iTrial) = ...
%         mean(phaseEGG(1,markers_downsampled_MRI(iTrial,3):markers_downsampled_MRI(iTrial,4)));
% end
% 
% % Control, check that averaging phase is working
% figure
% plot(angle(phaseEGG(markers_downsampled_MRI(1,3):end)),'-o')
% hold on
% plot(angle(EGG_phaseXVolume),'-or')
%
%% cut beginning & end of phase if necessary
% for n = 1:size(EGG_phaseXVolume.trial{1,1}, 1)
%     EGG_phaseXVolume.trial{1,1}(n,:) =EGG_phaseXVolume.trial{1,1}(n,(beginCut:endCut)); % cut begining and end of IRM acquisition
% end
% EGG_phaseXVolume.time{1,1} = EGG_phaseXVolume.time{1,1}(n,(beginCut:endCut));
