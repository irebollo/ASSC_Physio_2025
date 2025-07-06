%% general setup
clearvars
addpath('D:\ASSC28\_toolboxes\fieldtrip-20190329\');   % this old FieldTrip version has some nice features we want
ft_defaults;                                            % set defaults for FieldTrip

datapath = 'D:\ASSC28\_data\';                         % define where the data are
outpath = 'D:\ASSC28\_out\';                           % define where the output will go

ids = [1];                                              % demo sample (n = 1)

%% very basic, non-invasive preprocessing for M/EEG data
for isub = 1:numel(ids)
    
    % load raw data set
    load([datapath 'meg\rawmeg_vp' num2str(ids(isub)) '_rest_nasal.mat'],'raw');

    % resample everything to 300 Hz
    cfg             = [];               % start with an empty 'configuration' structure
    cfg.resamplefs  = 300;              % set 300 Hz as the target sampling rate
    rsraw = ft_resampledata(cfg,raw);   % run the command
    clear raw                           % clean up after yourself
     
    % apply DFT filter to remove line noise
    linefreq = 50;                                          % line frequency
    nharm = ceil((rsraw.fsample * 0.5) / linefreq) - 2;     % number of harmonics
    linefreqvec = linefreq:linefreq:(nharm+1)*linefreq;     % full vector of to-be-removed frequencies
    cfg                     = [];                           % empty structure
    cfg.channel             = {'MEG'};                      % apply to MEG channels only    
    cfg.dftfilter           = 'yes';                        % use DFT filter (which is the whole purpose of this call)
    cfg.demean              = 'yes';                        % execute demeaning
    cfg.dftfreq             = linefreqvec;                  % freqs of interest
    cfg.dftreplace          = 'neighbour';                  % replace by neighbouring frequencies
    cfg.dftbandwidth        = ones(1, nharm+1)*0.5;         % bandwidth of the line noise frequency
    cfg.dftneighbourwidth   = ones(1, nharm+1)*2;           % bandwidth of the neighbouring replacements
    megdft = ft_preprocessing(cfg,rsraw);                   % run the command
    clear rsraw
    
    % run ICA and mark eye-blink and cardiac components
    cfg             = [];                       % empty struct
    cfg.method      = 'runica';                 % run an ICA (instead of PCA, ...)
    cfg.runica.pca  = 20;                       % extract the first 20 components
    comps = ft_componentanalysis(cfg, megdft);  % run the command
    
    % visualise ICA results and select artefacts, if needed
    cfg             = [];               % empty structure
    cfg.layout      = 'CTF275.lay';     % layout to use: 275 sensors for CTF MEG scanners
    cfg.blocksize   = 8;                % look at 8 second windows
    cfg.channel     = 'all';            % show all components in one window
    ft_databrowser(cfg,comps);          % open the visualisation
    rejcomp = ft_icabrowser(cfg,comps); % open a second window in which we select the artificial components  
        
    %%% GETTING TO KNOW THE DATA %%%
    % Which components would you reject and why?             
    % Explore the data by comparing time series before and after running the ICA

    % remove selected components and save clean data
    cleanmeg = ft_rejectcomponent(struct('component',find(rejcomp)),comps,megdft);  % this command actually removes the components we marked above
    save([outpath 'cleanmeg_vp' num2str(ids(isub)) '_rest_nasal.mat'],'cleanmeg');  % save clean MEG data
end


cfg             = [];               % empty structure
cfg.layout      = 'CTF275.lay';     % layout to use: 275 sensors for CTF MEG scanners
cfg.blocksize   = 8;                % look at 8 second windows
cfg.channel     = 'all';            % show all components in one window
ft_databrowser(cfg,cleanmeg);          % open the visualisation
