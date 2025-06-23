function [FFT_EGG, figure_fft] = compute_FFT_EGG_peakchan(EGG_raw)

%%% Edited for Brain-Body Waves 2024 LB 
 
%% Create 200sec trials to run FFT on each trial 
% (resulting power spectra will be averaged across trials - giving a smooth power output) 

% cut data into trials by defining the length (in sec) of the data and the 
% overlap between segments (ratio)
cfg = [];
cfg.length                  = 200;
cfg.overlap                 = 0.75;
EGG_trials                  = ft_redefinetrial(cfg, EGG_raw);

%% Calculate spectrum

%Now we will run a Hann tapered fft on each of the trials, and the resulting power spectra will
% be averaged for us giving a smooth power output fft
cfg = [];
cfg.method = 'mtmfft'; % multitaper method with Fourier transform
cfg.taper = 'hanning'; % tapering window applied to the data before performing the Fourier transform is a Hanning window. 
cfg.output = 'pow'; % power spectrum
cfg.keeptrials      = 'no';
cfg.pad = 1000; % zero padding - seconds
cfg.foilim = [0 0.1]; % frequency range 0 - 6 cpm
FFT_EGG = ft_freqanalysis(cfg,EGG_trials);

%% find channel with highest peak in normogastric range

% Search for the largest peak within 2-4 cycles per min (normogastria)
% normogastria = 0.033-0.066 Hz (2-4 cpm or 30-15 secs)
indexNormogastria = find (FFT_EGG.freq >= 0.03333 & FFT_EGG.freq <= 0.06666); % normogastria window

%Automatically get the channel with highest peak in normogastria window
nChannels = length(EGG_raw.label);
maxPowerXChannel = zeros(nChannels,2); % column 1 max power, column 2 frequency location
for iChannel=1:nChannels
    maxPowerXChannel(iChannel,1) = max(FFT_EGG.powspctrm(iChannel,indexNormogastria));% max power within 0.033 to 0.066 hz window
    maxPowerLocation = FFT_EGG.powspctrm(iChannel,:)==maxPowerXChannel(iChannel,1);
    maxPowerXChannel(iChannel,2) = FFT_EGG.freq(find(maxPowerLocation)); % Freq of max peak
end

% find channel with maximum power (& frequency of that max power peak)
[highestPower, mostPowerfullChannel] = max(maxPowerXChannel(:,1));
mostPowerfullFrequency = maxPowerXChannel(mostPowerfullChannel,2);

%% Plot spectrum

% Create box for the normogastric win (of MaxPower for each channel)
NormogastriaPlot=zeros(nChannels,100);
for iChannel=1:nChannels
    NormogastriaPlot(iChannel,indexNormogastria)= maxPowerXChannel(iChannel,1);
end

% Create box for where gastric filter window will be (of MaxPower for each channel)
filterWidth= 0.0150;% Attention for visualization in the spectrum domain only, filter shape is different
filterPlot=zeros(nChannels,100);
for iChannel=1:nChannels
    indexFrequenciesFilter = find (FFT_EGG.freq >= maxPowerXChannel(iChannel,2)-filterWidth & FFT_EGG.freq <= maxPowerXChannel(iChannel,2)+filterWidth);
    filterPlot(iChannel,indexFrequenciesFilter)= maxPowerXChannel(iChannel,1);
end

% Which frequency range is going to appear in the plot
indexFrequenciesPloting = find (FFT_EGG.freq >= 0.0189 & FFT_EGG.freq <= 0.0698);

% plot FFT each channel
figure_fft= figure('visible','on');
for iChannel=1:nChannels
    subplot(round(nChannels/2),2,iChannel);
    plot(FFT_EGG.freq(indexFrequenciesPloting),FFT_EGG.powspctrm(iChannel,indexFrequenciesPloting),'-o','lineWidth',3);
    if iChannel == mostPowerfullChannel
        title(strcat(('Welch EGG n'),num2str(iChannel)),'fontweight','bold', 'fontsize',11);
    else
        title(strcat(('Welch EGG n'),num2str(iChannel)), 'fontsize',11);
    end
    hold on;
    plot (FFT_EGG.freq(indexFrequenciesPloting),filterPlot(iChannel,indexFrequenciesPloting),'r','lineWidth',3) % add EGG filter window in red
    plot (FFT_EGG.freq(indexFrequenciesPloting),NormogastriaPlot(iChannel,indexFrequenciesPloting),'k','lineWidth',3) % add normogastric win in black

    set(gca,'fontsize',11)
    xlim([0.0189 0.0698])

    set(gcf,'units','normalized','outerposition',[0 0 1 1])
    set(gcf, 'PaperPositionMode', 'auto');
end

%% manually select which channel to use (which channel has clearest & most powerful gastric peak)
message1=strcat(' Most powerfull channel is channel',32,num2str(mostPowerfullChannel));
fprintf(message1)

bestChannel = str2double(input('\nEnter the number of the channel you want to use:\n' ,'s'));
mostPowerfullFrequency = maxPowerXChannel(bestChannel,2);

message2=strcat(' Most powerfull frequency in this channel is ',32,num2str(mostPowerfullFrequency),'do you want to use that frequency?');
fprintf(message2)

useAutomaticFrequency = str2double(input('\n if yes enter 1,if Not put the number of the frequency you want to use:\n' ,'s'));
if useAutomaticFrequency ~=1
    mostPowerfullFrequency = useAutomaticFrequency;

    message3=strcat(' The  power is this channel is ',32,num2str(maxPowerXChannel(bestChannel,1))...
        ,32,' and the peak frequency is',32, num2str(mostPowerfullFrequency) );
    fprintf('\n')
    fprintf(message3)

    presstocontinue = str2double(input('\nPress any key to continue:\n' ,'s'));

end % not automatic channel selection

%% store important info
logEGGpreprocessing.mostPowerfullFrequency = mostPowerfullFrequency;
logEGGpreprocessing.selectedChannelPower = maxPowerXChannel(bestChannel,1);
logEGGpreprocessing.bestChannel = bestChannel;
logEGGpreprocessing.nChannels = nChannels;
logEGGpreprocessing.confidencechannel_best = str2double(input('\nHow confident you are this is the best channel ?(0 to 1): \n' ,'s'));
logEGGpreprocessing.confidencechannel_quality = str2double(input('\nHow you would rate the quality of the signal of that channel ?(0 to 1): \n' ,'s'));

FFT_EGG.logEGGpreprocessing = logEGGpreprocessing;
end


