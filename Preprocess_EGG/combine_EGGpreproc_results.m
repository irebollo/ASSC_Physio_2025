function [EGG_preproc] = combine_EGGpreproc_results(EGG_raw, EGG_phaseXVolume, MRI_TR, bestChannel)

cfg = [];
cfg.channel         = bestChannel;
EGG_raw             = ft_selectdata(cfg, EGG_raw);

disp('Resampling...to fMRI freq (TR)')
cfg = [];  %initialize configuration structure
cfg.detrend = 'no'; 
cfg.demean = 'yes';
cfg.resamplefs= 1/MRI_TR; 
EGG_downsampled_MRI = ft_resampledata(cfg,EGG_raw);

EGG_preproc = EGG_downsampled_MRI;
EGG_preproc.trial{1,1}(2:5,:) = EGG_phaseXVolume.trial{:}(:,:);
EGG_preproc.label = {'rawEGG', 'phase', 'amplitude', 'filteredEGG', 'phase_degrees'};

end