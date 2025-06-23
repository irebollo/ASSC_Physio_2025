function  error_csf_z=timeseries_csfSignal_obtainAndRegress(BOLD_filtered_zscored,cfgMain)

%{

This function retrieves the timeseries of a 3x3 cube of voxels located in
the CSF and regress it out the of all the brain voxels, 
the timeseries, the residuals of the regression and the betas coefficient are all stores in the timeseries
folder of each subject.
The rest of PHYSIENS analysis are performed on the residuals of this
regression, stored in the timeseries/data/ folder with the
csfRegressionResiduals_FB_S_SUBJECTNUMER name
The timeseries scripts have the suffix _regression to reflect. 

Output
Bold timeseries csf_regressed
Y:\Subjects\Subject13\Timeseries\MRItimeseries\csfRegressionResiduals_FB_S_13_kw3

IR commented the 12/09/2016
Checked 28/06/2017

%}

%% Input output filenames


%% Load files

gray_mask_filename = strcat(cfgMain.RootPath,'scripts',filesep,'files',filesep','gray_mask.hdr');

grayMatter = ft_read_mri(gray_mask_filename);
grayMatter = logical(grayMatter.anatomy);
grayMatter=grayMatter(:);

% tools_mni2vox([0 -46 -32]',3)
ventricles_coordinates = [27 23 7]; % Chosen manually from the 4th ventricle based on SPM apriori regions, mni coordinates are   0    -46   -32 
csf = zeros(53,63,46);

% Make a cube around the coordinates of csf
csf(ventricles_coordinates(1):ventricles_coordinates(1)+1,ventricles_coordinates(2):ventricles_coordinates(2)+1,ventricles_coordinates(3):ventricles_coordinates(3)+1) = 1;
% transform the 3d matrix to a vector
csf = csf(:);

insideBrain = tools_getIndexBrain('inside',cfgMain); % voxels inside brain



% REGRESSION

toBeExplained = BOLD_filtered_zscored'; % BOLD timeseries will be the variable to be explained out in the GLM
csf_timeseries = nanmean(BOLD_filtered_zscored(logical(csf),:)); % Average the timeseries in the csf compartment
explainingVariables = zscore([csf_timeseries]',[],1); % Variables used to explain the BOLD data
betas_csf = tools_EfficientGLM(toBeExplained,explainingVariables); % Obtain the betas indicating how much the predicting variable predicts the data
predictedBOLD = explainingVariables(:,1) *betas_csf(1,:); % What the BOLD timeseries should look like if CSF predicted at 100% accuracy the data
error_csf = toBeExplained - predictedBOLD; % The error is the portion of the data not predicted by the CSF signal
error_csf_z = zscore(error_csf,[],1);
error_csf_z = error_csf_z(:,insideBrain); 
betas_csf = betas_csf(:,insideBrain);


%% Plot regression to check it works



voxelCoordinates = sub2ind([53,63,46],11,30,37);
voxelCoordinates_inside = zeros(53*63*46,1);
voxelCoordinates_inside(voxelCoordinates)=1;
voxelCoordinates_inside = voxelCoordinates_inside(insideBrain);

if cfgMain.plotFigures == 0;
    SanityPlot = figure('visible','off');
else
    SanityPlot = figure('visible','on');
end

subj_idx= 99

hold on
subplot(2,1,1)
plot(csf_timeseries,'g','LineWidth',2)
title(['S',sprintf('%.2d',subj_idx),32,'CSF timeseries'],'fontsize',18)
grid on
subplot(2,1,2)
plot(toBeExplained(:,voxelCoordinates),'r-','LineWidth',4)
hold on
plot(error_csf_z(:,logical(voxelCoordinates_inside)),'b--','LineWidth',3)
grid on
legend ('Before regression','After regression')
title(['S',sprintf('%.2d',subj_idx),32,'EFFects of CSF regressionin right SS'],'fontsize',18)

   
set(gcf,'units','normalized','outerposition',[0 0 1 1])
set(gcf, 'PaperPositionMode', 'auto');





end