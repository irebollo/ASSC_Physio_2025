%% general setup
clearvars
datapath = 'D:\ASSC28\_data\';                                     % define where the data are
load([datapath 'resp\ASSC_resp_signal_preprocessing_raw.mat']);    % load the raw respiratory recording

%% normalizing and smoothing the data using a moving average (window length = 400ms)
% ! amount of smoothing can be adapted by choosing different window lengths
r = resp_raw_signal;
r(abs(zscore(r)) > 2.5) = NaN;                                                  % normalise time series and set outliers to NaN
r(isnan(r)) = interp1(find(~isnan(r)), r(~isnan(r)), find(isnan(r)), 'linear'); % linear interpolation of outlier segments
x = (r - nanmean(r))/nanstd(r);                                                 % (manually) normalise the interpolated time series
smresp = movmean(x,0.4*rsraw.fsample);                                          % light smoothing
zsmresp = zscore(smresp);  

% as a sanity check overlaying raw & smoothed signal
plot(timevec, data_resp_raw, 'LineWidth', 1.5, 'Color',[0.2 0.5 0.9])
hold on
plot(timevec,resp_smoothed_signal,'LineWidth', 2, 'Color',[0.2 0.2 0.2 0.8])
xlabel('time(sec)');
ylabel('au')
box off
set(gca,'TickDir','out', 'FontSize', 14, 'LineWidth', 1);
legend({'raw' 'smoothed'}, 'Location', 'southeast')
legend('boxoff')

%%% GETTING TO KNOW THE DATA %%%
% re-compute smoothing with different window lengths and compare the output signals

% double interpolation for phase extraction
    [~, peaks] = findpeaks(resp_smoothed_signal, 'MinPeakDistance',fs,'MinPeakProminence',0.1,'MinPeakHeight',0.1);
    troughs = [];
    for k = 2:numel(peaks)                  % start with peak #2
        tmp = zsmresp(peaks(k-1):peaks(k)); % get respiration course between peaks #1 and #2
        inds = peaks(k-1):peaks(k);         % get indices of the respiration course
        troughinds = inds(tmp == min(tmp)); % minimum of respiration between peaks is the trough
        troughs(k-1) = troughinds(1);       % for the rare but annoying case there is a peak longer than 1 sample
    end

    % compute phase vector
    phasevec = NaN(size(zsmresp));                          % make phase vector out of NaNs; length equals that of respiration vector
    phasevec(peaks) = 0;                                    % set phase at inspiration peaks to zero...
    phasevec(troughs) = pi;                                 % ... and troughs to pi
    for k = 1:numel(peaks)-1                                % loop over all peaks
        tmp = phasevec(peaks(k)+1:troughs(k));              % phase angles between peak #1  and trough #1
        tmp2 = phasevec(troughs(k)+1:peaks(k+1));           % phase angles between trough #1 and peak #2
        sub = linspace(0+pi/numel(tmp),pi,numel(tmp));      % linear interpolation peak2trough
        sub2 = linspace(-pi+pi/numel(tmp),0,numel(tmp2));   % same for trough2peak
        phasevec(peaks(k)+1:troughs(k)) = sub;              % substitute with phase angles (peak2trough)
        phasevec(troughs(k)+1:peaks(k+1)) = sub2;           % same for trough2peak
    end

    % add everything to the plot
    plot(zsmresp,'k');                          % plot respiration in black
    hold on 
    scatter(peaks,zsmresp(peaks),10,'r');       % overlay peaks in red...
    scatter(troughs,zsmresp(troughs),10,'b');   % ... and troughs in blue
    ylabel('resp');                             % add vertical label on the left
    xlim([0 length(zsmresp)]);                  % limit to length of recording
    ylim([-4 4]);                               % limits for z score for resp signal
    yticks([-4 0 4]);                           % set ticks for left side
    yyaxis right                                % right side now
    plot(phasevec,'Color',[0 0 0 0.2]);         % add phase vector
    ylabel('resp phase');                       % add vertical label on the right
    ylim([-pi pi]);                             % limits for resp phase
    yticks([-pi 0 pi]);                         % set ticks for right side
    yticklabels({'-pi' '0' 'pi'});              % add tick labels
    xlabel('time (samples)');                   % add x axis label on bottom plot
    hold off

%%% GETTING TO KNOW THE DATA %%%
% compare the phase vector you just computed to what a simple Hilbert transform would give you
