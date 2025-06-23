function BOLD_filtered_zscored = timeseries_preprocessBOLD(BOLDtimeseries,cfgMain)

%{

Takes the output of timeseries_prepare_import2matlab that is stored in 
each subject timeseries\data folder 
Y:\Subjects\Subject13\Timeseries\MRItimeseries\fMRItimeseries_S13_kw3.mat

and preprocess it so it can be further analyzed
Preprocessing steps, implemented in fieldtrip, includes remove polinomial trends of second degree,
bandpass filter between 0.01 and 0.1 Hz and standarize to Z units
the output is stored in the tiemseries folder of each subject
Y:\Subjects\Subject13\Timeseries\MRItimeseries\BOLDFULLBAND_S13_kw3


% IR commented on the 12/09/2016
% recommented 28/06/2017

%}



BOLDtimeseries.BOLDtsNOPOLI = ft_preproc_polyremoval (BOLDtimeseries.trialVector', 2);
BOLDtimeseries.BOLD_filtered = ft_preproc_bandpassfilter (BOLDtimeseries.BOLDtsNOPOLI,0.5,[0.01 0.1]);
BOLD_filtered_zscored = ft_preproc_standardize (BOLDtimeseries.BOLD_filtered);



%% Sanity check


    
    if cfgMain.plotFigures == 0;
        SanityPlot = figure('visible','off');
    else
        SanityPlot = figure('visible','on');
    end
    
    voxelCoordinates = sub2ind([53,63,46],11,30,37);
% voxelCoordinates = sub2ind([79,95,79],9,45,53); % voxel in somatomotor cortex
    subj_idx=99

    time=BOLDtimeseries.time;
    subplot(4,1,1) % 
    plot(time,BOLDtimeseries.trialVector(:,voxelCoordinates),'LineWidth',4)
    xlabel('seconds')
    title(['S',sprintf('%.2d',subj_idx),32,'raw'],'fontsize',18)
    subplot(4,1,2)
    plot(time,BOLDtimeseries.BOLDtsNOPOLI(voxelCoordinates,:),'LineWidth',4)
        xlabel('seconds')

    title(['S',sprintf('%.2d',subj_idx),32,'Polinomial removal'],'fontsize',18)
    subplot(4,1,3)
    plot(time,BOLDtimeseries.BOLD_filtered(voxelCoordinates,:),'LineWidth',4)
        xlabel('seconds')

    title(['S',sprintf('%.2d',subj_idx),32,'Bandpass filtered between 0.01 and 0.1'],'fontsize',18)
    subplot(4,1,4)
    plot(time,BOLD_filtered_zscored(voxelCoordinates,:),'LineWidth',4)
    title('zscored','fontsize',18)
        xlabel('seconds')

    
    set(gcf,'units','normalized','outerposition',[0 0 1 1])
    set(gcf, 'PaperPositionMode', 'auto');
    
  


end