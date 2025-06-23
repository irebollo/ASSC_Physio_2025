%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% SCRIPT_EGG_main_ASSC2025
%%% 
%%% Edited for ASSC 2025 IR
%%% Edited for Brain-Body Waves 2024 LB 
%%%
%%% This script calls the functions needed for EGG preprocessing, artifact 
%%% detection and EGG analysis.
%%% When using this function in any published study, please cite: Wolpert, 
%%% N., Rebollo, I., Tallon-Baudry, C. (2020). Electrogastrography for 
%%% psychophysiological research: practical considerations, analysis pipeline 
%%% and normative data in a large sample. Psychophysiology (in press)
%%%
%%% The scripts and functions were written in Matlab version R2017b.
%%%
%%% These functions make use of the fieldtrip toolbox, version 20170315
%%% (see http://www.fieldtriptoolbox.org/).
%%% Reference:
%%% Robert Oostenveld, Pascal Fries, Eric Maris, and Jan-Mathijs Schoffelen. 
%%% FieldTrip: Open Source Software for Advanced Analysis of MEG, EEG, and 
%%% Invasive Electrophysiological Data. Computational Intelligence and 
%%% Neuroscience, vol. 2011, Article ID 156869, 9 pages, 2011. 
%%% doi:10.1155/2011/156869.
%%% 
%%% The code comes with example of EGG datasets (EGG_raw_example1/2/3.mat).
%%% The files contains 7 channels of EGG recorded for approximately 12 
%%% minutes using a Biosemi amplifier (sampling rate: 1kHz).
%%%
%%% Copyright (C) 2019, Laboratoire de Neurosciences Cognitives, Nicolai 
%%% Wolpert, Ignacio Rebello & Catherine Tallon-Baudry
%%% 
%%% DISCLAIMER:
%%% This code is provided without explicit or implicit guarantee, and 
%%% without any form of technical support. The code is not intended to be 
%%% used for clinical purposes. The functions are free to use and can be
%%% redistributed, modified and adapted, under the terms of the CC BY-NC-SA
%%% version of creative commons license (see
%%% <https://creativecommons.org/licenses/>).
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% 1) Initialization
clear all; close all;
clc

% specify location of scripts and data files
script_path = '/mnt/data/Work/Stomach_Brain_ASSC_2025/EGG_Scripts-master/';
file_name   = 'EGGdata/EGG_raw_example1.mat';

% add fieldtrip toolbox to path
fieldtrip_path = '/mnt/data/Work/PresentationWavesPractical/fieldtrip/';
addpath(fieldtrip_path);
ft_defaults;
% add EGG functions to path 
addpath(genpath(script_path));

% Load EGG raw data.
% A dataset in fieldtrip format (e.g. output of 'ft_preprocessing') is expected here. 
% For an overview of fieldtrip-compatible dataformats see:
% http://www.fieldtriptoolbox.org/faq/dataformat/
% 
% As a reminder, EGG data should be recorded without highpass-filter (DC recording)
% and referenced in an appropriate manner. 
% For more information on reading, filtering and re-referencing with fieldtrip, 
% see: http://www.fieldtriptoolbox.org/tutorial/continuous/
load(strcat([script_path filesep file_name]));

%cfg = [];
%cfg.dataset = '/Users/au704655/Documents/StomachBrain/WAVES2024/EGG_Scripts-master/EGGdata/EGG_raw_example4.set';
%EGG_raw = ft_preprocessing(cfg)

%% downsample data to 10 Hz

disp('Resampling...')
cfg = [];  %initialize configuration structure
cfg.detrend = 'no'; 
cfg.demean = 'yes';
cfg.resamplefs= 10; % 10 Hz
EGG_downsampled = ft_resampledata(cfg,EGG_raw); % This procedure also lowpass filter the data at half the new sr

%% plot raw EGG time-series

[figure_EGGtimeseries] = plot_EGGtimeseries(EGG_downsampled);

%% Power spectrum, channel selection
% Compute the EGG power spectrum for all channels and select the channel
% with the strongest peak between 0.033 and 0.067 Hz.

%[FFT_EGG,figure_fft] = compute_FFT_EGG(EGG_downsampled);
[FFT_EGG, figure_fft] = compute_FFT_EGG_peakchan(EGG_downsampled); % manual channel selection

%% Filtering around gastric peak frequency
% Filter the raw EGG from the selected channel around the peak dominant
% frequency to extract the gastric rhythm using a finite impulse response
% filter (Matlab FIR2), with a banwith of +/- 0.015 Hz of the peak
% frequency.

%EGG_filtered = compute_filter_EGG(EGG_downsampled, FFT_EGG); % also computes phase if interested EGG processing without fMRI
EGG_filtered = compute_filter_EGG_peakchan(EGG_downsampled, FFT_EGG);

%% Compute Average EGG phase per MRI volume

MRI_TR = 1.4;
EGG_phaseXVolume = compute_phaseXvolume(EGG_filtered, MRI_TR);

%% Combine EGG preproc results 

EGG_preproc = combine_EGGpreproc_results(EGG_downsampled, EGG_phaseXVolume, MRI_TR, FFT_EGG.logEGGpreprocessing.bestChannel);

%% Visual inspection

[figure_EGG_preproc] = plots_EGG_visualinspection(EGG_downsampled, FFT_EGG, EGG_phaseXVolume, FFT_EGG.logEGGpreprocessing.bestChannel);

%% Show standard deviation of cycle duration

stds_cycle_durations = compute_std_cycle_duration(EGG_preproc);

%% Show proportion of normogastria

show_prop_normogastria(EGG_preproc)

%% Artifact detection

close all;
[art_def, figure_EGG_artifacts] = detect_EGG_artifacts(EGG_preproc);
