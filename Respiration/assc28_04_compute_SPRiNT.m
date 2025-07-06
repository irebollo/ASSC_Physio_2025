%% general setup
clearvars
addpath('D:\ASSC28\_toolboxes\\SPRiNT\');  % add the SPRiNT toolbox

datapath = 'D:\ASSC28\_data\';     % define where the data are
outpath = 'D:\ASSC28\_out\';       % define where the output will go

ids = [1];    % demo sample (n = 1)

%% using the SPRiNT toolbox, we will compute time-resolved 1/f slope and power spectra
for isub = 1:numel(ids)
    load([outpath 'cleanmeg_vp' num2str(ids(isub)) '_rest_nasal.mat'],'cleanmeg');      % clean MEG data

    % set up SPRiNT parameters
    % input
    F = cleanmeg.trial{1}(1,:);         % input time series - single channel for demo purposes

    % STFT opts
    opt.sfreq = 300;                    % input sampling rate
    opt.WinLength = 1;                  % 1; STFT window length
    opt.WinOverlap = 75;                % overlap between sliding windows (in %)
    opt.WinAverage = 1;                 % 5; number of sliding windows averaged by time point

    % specparam opts
    opt.freq_range          = [1 40];
    opt.peak_width_limits   = [0.5 6];  % default from the paper
    opt.max_peaks           = 3;
    opt.min_peak_height     = 6 / 10;   % convert from dB to B
    opt.aperiodic_mode      = 'fixed';  % alternative: knee
    opt.peak_threshold      = 2.0;      % 2 std dev: parameter for interface simplification

    % Matlab-only options
    opt.peak_type           = 'gaussian';   % alternative: cauchy
    opt.proximity_threshold = 2;
    opt.guess_weight        = 'none';
    opt.thresh_after        = true;         % threshold after fitting, always selected for Matlab
    % (mirrors the Python FOOOF closest by removing peaks
    % that do not satisfy a user's predetermined conditions)
    % only used in the absence of the
    if license('test','optimization_toolbox') % check for optimization toolbox
        opt.hOT = 1;
        disp('Using constrained optimization, Guess Weight ignored.')
    else
        opt.hOT = 0;
        disp('Using unconstrained optimization, with Guess Weights.')
    end
    opt.rmoutliers          = 'yes';
    opt.maxfreq             = 2.5;
    opt.maxtime             = 6;
    opt.minnear             = 3;

    Freqs = 0:1/opt.WinLength:opt.sfreq/2;
    channel = struct('data',[],'peaks',[],'aperiodics',[],'stats',[]);
    % Compute short-time Fourier transform
    [TF, ts] = SPRiNT_stft(F,opt);
    outputStruct = struct('opts',opt,'freqs',Freqs,'channel',channel);
    % Parameterize STFTs
    s_data = SPRiNT_specparam_matlab(TF,outputStruct.freqs,outputStruct.opts,ts);

    % save the output
    save([outpath 'sprint_vp' num2str(ids(isub)) '.mat'],'s_data');
end

%%% GETTING TO KNOW THE DATA %%%
% Explore the vast amount of data stored in the s_data structure.
% What kinds of information can you gather from there?

%% one example of how to relate time-resolved 1/f slope to respiration phase
nbin = 60;                      % number of phase bins to use
phw = 2*pi/10;                  % width of each phase bin 
pb = linspace(-pi,pi,nbin);     % borders of the phase bins

binslopes = zeros(numel(ids),nbin); % start with a zero structure

for isub = 1:numel(ids)
    load([outpath 'sprint_vp' num2str(ids(isub)) '.mat'],'s_data');                     % SPRiNT output
    load([outpath 'phasevec_vp' num2str(ids(isub)) '_rest_nasal.mat'],'phasevec');      % resp phase

    % get respiration phase at SPRiNT bin centres (ie time points around which 1/f slope was computed
    binctrs = [s_data.SPRiNT.channel.data.time];    % time points 
    binresp = phasevec(binctrs*300);                % corresponding respiration phase

    % work through respiration and compute average 1/f slope
    statslope = [s_data.SPRiNT.channel.aperiodics.exponent];    % for now, we look at 1/f slope only

    for ib = 1:nbin                                                     % loop through all phase bins
        inds = find((binresp>pb(ib)-phw) & (binresp<pb(ib)+phw) | ...   % find all 'trials' (= slope estimates) falling within this phase bin
            (binresp-2*pi>pb(ib)-phw) & (binresp-2*pi<pb(ib)+phw)| ...
            (binresp+2*pi>pb(ib)-phw) & (binresp+2*pi<pb(ib)+phw));
        binslopes(isub,ib) = mean(statslope(inds));                     % average all estimates within this phase bin
    end
end

%%% GETTING TO KNOW THE DATA %%%
% Plot data within the binslopes structure; also normalise the data across
% the respiration phase dimension to remove baseline effects.
% Re-run the analysis with changed parameters of number of bins (nbins) and
% width of the phase bins. Try making it as simple as inspiration vs
% expiration. 
