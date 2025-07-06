%% general setup
clearvars

datapath = 'D:\ASSC28\_data\';     % define where the data is
outpath = 'D:\ASSC28\_out\';       % define where the output will go

% full sample for demo
ids = [2549 3051 3066 3048 3033 3050 3067 3064 3062 3063 2521 2550 2448 2314 2168 ...
       2524 2480 2299 2749 1566 2214 2348 1822 3090 3047 3081 3053 2330 3078 3105]; 

t = readtable([datapath 'behav\fullbehavtable.txt']);     % load one big table with all the results 

%% start with some very basic behavioural measures
% overall hit rate
for isub = 1:numel(ids)                                             % loop over participants
    tmp = t(t.ID == ids(isub),:);                                   % find data for this participant
    fullhr(isub) = length(find(tmp.corr == 1))/height(tmp);         % compute hit rate
end
m = mean(fullhr)
sd = std(fullhr)

% hit rate across conditions
allhr = zeros(numel(ids),4);                                                % empty struct, now for the 4 conditions separately
for isub = 1:numel(ids)                                                     % loop over participants
    tmp = t(t.ID == ids(isub),:);                                           % find data for this participant
    for icond = 1:4                                                         % loop over conditions
        tmp2 = tmp(tmp.cond == icond,:);                                    % take only this condition
        allhr(isub,icond) = length(find(tmp2.corr == 1))/height(tmp2);      % condition-specific hit rate
    end
end

[p,tbl,stats] = friedman(allhr);       % significant differences between conditions?
[results,m,h,gnames] = multcompare(stats)   % where are the differences?

%%% GETTING TO KNOW THE DATA %%%
% Explore the results using the interactive last two lines of code.
% What do the stats tell you?

%% hit rate per condition over respiration
t.corr(t.corr == 9) = 0;    % this counds false alarms as incorrect

nsub = length(ids);             % number of participants
nbin = 60;                      % number of resp phase bins
phw = 2*pi/10;                  % width of the phase bins
pb = linspace(-pi,pi,nbin);     % borders of the phase bins

binhr = zeros(4,numel(ids),nbin);           % empty struct
for isub = 1:nsub                           % loop over participants
    tmp = t(t.ID == ids(isub),:);           % this participant only
    for icond = 1:4                         % loop over conditions
        tmp2 = tmp(tmp.cond == icond,:);    % select condition
        for ib = 1:nbin  	                                                                % loop over bins
            inds = find((tmp2.stimphase>pb(ib)-phw) & (tmp2.stimphase<pb(ib)+phw) | ...     % find trials falling into this bin
                (tmp2.stimphase-2*pi>pb(ib)-phw) & (tmp2.stimphase-2*pi<pb(ib)+phw)| ...
                (tmp2.stimphase+2*pi>pb(ib)-phw) & (tmp2.stimphase+2*pi<pb(ib)+phw));
            tmp3 = tmp2(inds,:);                                                            % take data from all trials within this bin
            binhr(icond,isub,ib) = length(find(tmp3.corr == 1))/height(tmp3);               % compute bin-wise hit rate
        end
    end
end

% plot
conds = {'C+T+','C-T+','C+T-','C-T-'};
cols = {'#048789','#503D2E','#D44D27','#E2A72E'};
mzbinhr = squeeze(median(zscore(binhr,[],3),2));
figure;
yline(0,'k--');
hold all
for icond = [1 4]
    plot(mzbinhr(icond,:),'Color',[hex2rgb(cols{icond}) 0.5], 'LineWidth',1);
end
xlabel('phase bin');
ylabel('z(hit rate)');
legend('',conds{1},conds{4})
ylim([-1 1]);
title('HR ~ resp phase per condition');

%%% GETTING TO KNOW THE DATA %%%
% How would you summarise the main result in one sentence? 
% Re-run the analysis with changed parameters of number of bins (nbins) and
% width of the phase bins. Try making it as simple as inspiration vs
% expiration. 

%% finally, let's run a linear mixed-effects model on hit rate
pb = linspace(-pi,pi,60);

% set up a smaller table on which we'll run the LMEM
lmmtbl  = [];   % empty LMEM table
subj    = [];   % participant vector
cond    = [];   % condition
ph1     = [];   % phase vector 1 - cosine
ph2     = [];   % phase vector 2 - sine
hr      = [];   % hit rate

for k = 1:numel(ids)                                % loop over participants
    for k2 = 1:4                                    % loop over conditions
        subj    = [subj k*ones(1,60)];              % fill in the participant ID
        cond    = [cond k2*ones(1,60)];             % condition
        ph1     = [ph1 cos(pb)];                    % cosine
        ph2     = [ph2 sin(pb)];                    % sine
        hr      = [hr squeeze(binhr(k2,k,:))'];     % hit rate
    end
end

varnames = {'subj','cond','cos' 'sin' 'hitrate'};                           % variable names
lmmtbl = table(subj', cond', ph1', ph2', hr', 'VariableNames',varnames);    % construct the table

model1 = 'hitrate ~ cond + cos + sin + (1|subj)';                                           % define the model
lme1 = fitlme(lmmtbl,model1);                                                               % fit the model
truebeta = sqrt(lme1.Coefficients(3,2).Estimate^2 + lme1.Coefficients(4,2).Estimate^2 );    % get the vector norm of sine and cosine betas

% null distribution
for k3 = 1:1000                     % loop over randomisation iterations
    disp(['n = ' num2str(k3)]);     % how much longer?
    
    lmmtbl  = [];
    subj    = [];
    cond    = [];
    ph1     = [];
    ph2     = [];
    hr      = [];
    
    for k = 1:numel(ids)
        for k2 = 1:4
            subj    = [subj k*ones(1,60)];
            cond    = [cond k2*ones(1,60)];
            ph1     = [ph1 cos(pb)];
            ph2     = [ph2 sin(pb)];
            xx      = squeeze(binhr(k2,k,:));
            hr      = [hr xx(randperm(size(xx,1)))'];   % random permutation of hit rates across bins
        end
    end
    
    varnames = {'subj','cond','cos' 'sin' 'hitrate'};
    lmmtbl = table(subj', cond', ph1', ph2', hr', 'VariableNames',varnames);
    
    model1 = 'hitrate ~ cond + cos + sin + (1|subj)';
    lme1 = fitlme(lmmtbl,model1);
    nullbetas(k3) = sqrt(lme1.Coefficients(3,2).Estimate^2 + lme1.Coefficients(4,2).Estimate^2 );   % save the vector norm beta each time
end

% get percentile of empirical vector norms
[perm_ecdf(:, 1), perm_ecdf(:, 2)] = ecdf(nullbetas);       % empirical density function of the 'null betas'
last_idx = find(perm_ecdf(:, 2) < truebeta, 1, 'last');     % locate the true beta relative to that density function
p = 1-perm_ecdf(last_idx)                                   % compute p value 

% plot
h = histfit(nullbetas,100,'kernel');
h(1).FaceColor = [255 255 255]./255;
h(1).EdgeColor = [0 0 0];
h(2).Color = [[220 53 34]./255 0.75];
xline(truebeta, '-b', {['p = ' num2str(p)]});
axis('square')
xlabel('beta');
ylabel('# of obs');

%%% GETTING TO KNOW THE DATA %%%
% Think of other model configurations that could be interesting.
% How do the results change with different model definitions?
% If you're unfamiliar with LMEMs, try 'help fitlme'
