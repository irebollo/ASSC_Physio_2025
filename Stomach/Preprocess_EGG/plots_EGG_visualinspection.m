function figure_EGG_preproc = plots_EGG_visualinspection(EGG_downsampled, FFT_EGG, EGG_phaseXVolume, bestChannel)

%%% Edited for Brain-Body Waves 2024 LB 

% This function plots the EGG raw signal, filtered signal, phase and amplitude
%
% Inputs
%     EGG_filtered          filtered EGG signal (output from 'compute_filter_EGG.m')
% 
% Outputs
%     figure_EGG_filtered figure showing EGG raw signal, filtered signal,
%                         phase and amplitude
%
% When using this function in any published study, please cite: 
% Wolpert, N, Rebollo, I, Tallon‐Baudry, C. Electrogastrography for psychophysiological 
% research: Practical considerations, analysis pipeline, and normative data in a large 
% sample. Psychophysiology. 2020; 57:e13599. https://doi.org/10.1111/psyp.13599 
%
% This function was written in Matlab version R2017b.
%
% This function make use of the fieldtrip toolbox, version 20170315
% (see http://www.fieldtriptoolbox.org/).
% Reference:
% Robert Oostenveld, Pascal Fries, Eric Maris, and Jan-Mathijs Schoffelen. 
% FieldTrip: Open Source Software for Advanced Analysis of MEG, EEG, and 
% Invasive Electrophysiological Data. Computational Intelligence and 
% Neuroscience, vol. 2011, Article ID 156869, 9 pages, 2011. 
% doi:10.1155/2011/156869.
%
% Copyright (C) 2019, Laboratoire de Neurosciences Cognitives, Nicolai 
% Wolpert, Ignacio Rebello & Catherine Tallon-Baudry
% Email: nicolaiwolpert@gmail.com
% 
% DISCLAIMER:
% This code is provided without explicit or implicit guarantee, and without 
% any form of technical support. The code is not intended to be used for 
% clinical purposes. The functions are free to use and can be 
% redistributed, modified and adapted, under the terms of the CC BY-NC-SA
% version of creative commons license (see
% <https://creativecommons.org/licenses/>).

% plot FFT
figure_EGG_preproc = figure('units','normalized','outerposition',[0 0 1 1]); set(gcf,'color','w');
subplot(5,1,1);
NormoGastricindexFrequencies = find (FFT_EGG.freq >= 0.035 & FFT_EGG.freq <= 0.065);
plot(FFT_EGG.freq(NormoGastricindexFrequencies),FFT_EGG.powspctrm(bestChannel,NormoGastricindexFrequencies),'-o','lineWidth',3);
ax = gca;
ax.XLim            = [FFT_EGG.freq(NormoGastricindexFrequencies(1)) FFT_EGG.freq(NormoGastricindexFrequencies(end))];
ax.FontSize        = 16;
ax.YLabel.String   = 'Power';
ax.YLabel.FontSize = 17;
title(sprintf('Power Spectrum: peak freq = %.3f Hz', FFT_EGG.logEGGpreprocessing.mostPowerfullFrequency), 'FontSize', 25);

% Plot EGG raw signal
subplot(5,1,2); 
plot(EGG_downsampled.time{1}, EGG_downsampled.trial{1}(bestChannel,:));
ax = gca;
ax.XLim            = [EGG_downsampled.time{1}(1) EGG_downsampled.time{1}(end)];
ax.FontSize        = 16;
ax.YLabel.String   = 'µV';
ax.YLabel.FontSize = 17;
title('Raw signal', 'FontSize', 25);

% Plot filtered signal
subplot(5,1,3);
plot(EGG_phaseXVolume.time{1}, EGG_phaseXVolume.trial{1}(4, :), 'g', 'LineWidth', 1);
ax = gca;
ax.XLim            = [EGG_phaseXVolume.time{1}(1) EGG_phaseXVolume.time{1}(end)];
ax.FontSize        = 16;
ax.YLabel.String   = 'µV';
ax.YLabel.FontSize = 17;
title('Filtered signal', 'FontSize', 25);

% Plot signal phase in rads
subplot(5,1,4);
plot(EGG_phaseXVolume.time{1}, EGG_phaseXVolume.trial{1}(1, :), 'g', 'LineWidth', 1);
ax = gca;
ax.XLim            = [EGG_phaseXVolume.time{1}(1) EGG_phaseXVolume.time{1}(end)];
ax.FontSize        = 16;
ax.YLabel.String   = 'Phase (rad)';
ax.YLabel.FontSize = 17;
title('Phase', 'FontSize', 25);

% Plot signal amplitude
subplot(5,1,5);
plot(EGG_phaseXVolume.time{1}, EGG_phaseXVolume.trial{1}(2, :), 'r', 'LineWidth', 1);
ax = gca;
ax.XLim            = [EGG_phaseXVolume.time{1}(1) EGG_phaseXVolume.time{1}(end)];
ax.FontSize        = 16;
ax.XLabel.String   = 'Time (seconds)';
ax.XLabel.FontSize = 17;
ax.YLabel.String   = 'µV';
ax.YLabel.FontSize = 17;
title('Amplitude', 'FontSize', 25);

end
