%% general setup
clearvars
addpath('D:\ASSC28\_code\');       % add the code path
datapath = 'D:\ASSC28\_data\';     % define where the data are
outpath = 'D:\ASSC28\_out\';       % define where the output will go

ids = [1];                          % demo sample (n = 1)
    
%% first, we need vectors of respiratory phase
for isub = 1:numel(ids)

    load([datapath 'resp\zresp_vp' num2str(ids(isub)) '_rest_nasal.mat'],'zresp'); % z-scored respiration signal

    % find peaks and troughs in the respiration time series
    [~, peaks] = findpeaks(zresp, 'MinPeakProminence',1); % the last argument uses a criterion of z = 1 for peak prominence
    troughs = [];
    for k = 2:numel(peaks)                  % start with peak #2
        tmp = zresp(peaks(k-1):peaks(k));   % get respiration course between peaks #1 and #2
        inds = peaks(k-1):peaks(k);         % get indices of the respiration course
        troughinds = inds(tmp == min(tmp)); % minimum of respiration between peaks is the trough
        troughs(k-1) = troughinds(1);       % for the rare but annoying case there is a peak longer than 1 sample
    end

    phasevec = NaN(size(zresp));            % initiate phase vector with NaNs 
    phasevec(peaks) = 0;                    % inspiration maximum = phase zero
    phasevec(troughs) = pi;                 % inspiration minimum = phase +/- pi
    for k = 1:numel(peaks)-1
        tmp = phasevec(peaks(k)+1:troughs(k));              % phase angles between peak #1  and trough #1
        tmp2 = phasevec(troughs(k)+1:peaks(k+1));           % phase angles between trough #1 and peak #2
        sub = linspace(0+pi/numel(tmp),pi,numel(tmp));      % linear interpolation peak2trough
        sub2 = linspace(-pi+pi/numel(tmp),0,numel(tmp2));   % same for trough2peak
        phasevec(peaks(k)+1:troughs(k)) = sub;              % substitute with phase angles (peak2trough)
        phasevec(troughs(k)+1:peaks(k+1)) = sub2;           % same for trough2peak
    end
    
    save([outpath 'phasevec_vp' num2str(ids(isub)) '_rest_nasal.mat'],'phasevec');
end

%%% GETTING TO KNOW THE DATA %%%
% Re-run the peak detection step with different parameter settings. 
% To what extent do the results change, and how?

%% compute the modulation index
for isub = 1:numel(ids)
    load([outpath 'phasevec_vp' num2str(ids(isub)) '_rest_nasal.mat'],'phasevec');      % resp phase
    load([outpath 'cleanmeg_vp' num2str(ids(isub)) '_rest_nasal.mat'],'cleanmeg');      % clean MEG data
    load([datapath 'hm\headmovtc_vp' num2str(ids(isub)) '_rest_nasal.mat'],'headmov');  % head motion time series

    % set up the head motion GLM
    hm = headmov.trial{1};
    hm = detrend(hm',3)';                   % remove third order polynomial
    hm = [hm; hm.^2; hm.^3];                % get derivatives and nonlinear regressors
    hm = [hm; [zeros(18,1) diff(hm,1,2)]];  % add derivatives
    hm = zscore(hm,0,2);                    % normalise along the second dimension
    hm = movmean(hm',10)';                  % very slight temporal smoothing (moving window average across 10 samples)

    gfp = zeros(60,cleanmeg.sampleinfo(2)); % set up empty GFP array - sixty is the number frequencies in the filter bank below
    fb = cwtfilterbank('SignalLength',cleanmeg.sampleinfo(2),'SamplingFrequency',cleanmeg.fsample,'FrequencyLimits',[2 120],'WaveletParameters',[3,20]);

    for ichan = 1:5                                                 % five channels for demo purposes
        [Envc,f] = wt(fb,cleanmeg.trial{1}(ichan,:));               % compute envelope per frequency (output is complex)
        y = abs(Envc);                                              % make complex into natural numbers
        y = flip(y);                                                % sort frequencies in ascending order (they are reversed in the filter bank)
        
        % regress out head movement
        y2 = y;                         % make a copy of the envelope (nfreqs x nsamples)
        for k = 1:length(f)             % loop across frequencies
            me = mean(y(k,:));          % get overall mean 
            offset = y(k,:)-me;         % get mean-centred envelope course
            b = hm'\offset';            % get 'weights' of head motion at each sample
            model = hm'*b;              % multiply head motion time course with its weights
            y2(k,:) = offset-model'+me; % remove head motion-related part of mean centred envelopes
        end
        
        gfp = gfp + abs(y2);            % simply sum up across channels
    end

    mgfp = movmean(gfp', 300)';         % slight temporal smoothing (moving average with a 1s window)
    r(isub,:) = mod_indexR2(phasevec(~isnan(phasevec)), mgfp(:,~isnan(phasevec))'); % compute modulation index

    % compute phase-triggered averages (= TFR over respiratory phase)
    zc = find(phasevec == 0);       % find zero crossings within the phase vector
    pta = zeros(length(f),2000);    % set up empty struct for PTA
    for k = 1:length(zc)        
        if ((zc(k)-999>0) &&  (zc(k)+1000<length(gfp)))     % get 2000 samples (~6 seconds) around each peak
            pta = pta + gfp(:,[zc(k)-999:zc(k)+1000]);      % get GFP at those samples
        end
    end
    normpta(isub,:,:) = normalize(pta,2); % normalise within each frequency
end

%%% GETTING TO KNOW THE DATA %%%
% Explore the outputs through some plots. What do observe in the PTA? 
% Also re-run the GFP and PTA analyses without correcting for head motion.
% To what extent do the results change?













