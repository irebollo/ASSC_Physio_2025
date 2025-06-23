%% Paths

addpath('/mnt/data/Work/Stomach_Brain_ASSC_2025/Coupling/EGG_Scripts-master')
addpath('/mnt/data/Work/PresentationWavesPractical/fieldtrip');ft_defaults

RootPath='/mnt/data/Work/Stomach_Brain_ASSC_2025/Coupling/';
addpath(genpath([RootPath,filesep,'scripts']));

datapath=[RootPath filesep 'data' filesep];

cfgMain=global_getcfgmain;
cfgMain.plotFigures=1;
cfgMain.RootPath=RootPath;
%% Load EGG data



% Load BOLD data in fieldtrip format
load([datapath 'fMRItimeseries_kw3.mat']);

BOLD_filtered_zscored=timeseries_preprocessBOLD(BOLDtimeseries,cfgMain); % preprocess the timeseries (remove polinomial and filter)
error_csf_z=timeseries_csfSignal_obtainAndRegress(BOLD_filtered_zscored,cfgMain);clear BOLD_filtered_zscored % obtain the wm and csf signals from the BOLD preprocessed timeseries and stores them
phaseMRI=timeseries_preparePhases_Regression(error_csf_z,cfgMain);clear error_csf_z % Filter and hilbert transform residuals of csf regression
[PLVMAP]=timeseries_mapPLV_Regression(phaseMRI,cfgMain);

[SURRPLV]=timeseries_medianRotation_Regression(phaseMRI,cfgMain);


%% Get coupling stength
medianRotationFilename=
couplingStrength = PLVMAP -SURRPLV;
tools_writeMri(couplingStrength,medianRotationFilename)

[insideBrain] = tools_getIndexBrain('inside',cfgMain);

subj_idx=99;
figure
nhist(couplingStrength(insideBrain))
xlabel('PLV')
title(['S',sprintf('%.2d',subj_idx),32,'Coupling strength across bain. Mean:' num2str(nanmean(couplingStrength(insideBrain))) ' rSS voxel:' 32 num2str(medianPLV(insideBrain(ind_voxelCoordinates_inside)))],'fontsize',18)
