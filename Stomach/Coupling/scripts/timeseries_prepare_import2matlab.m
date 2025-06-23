function timeseries_prepare_import2matlab(subj_idx,cfgMain)

%{
Construct fMRI timeseries (time,voxels) from swaf images located at each subject folder
and saves them as a structure in the HDD.

Input files
hdr files from each volume
Y:\Subjects\Subject13\fMRI\acquisition1\RestingState\s3wafPHYSIENS_Sujet13-0002-00001-000001-01.hdr

Output 
concatenated MRI timeseries
Y:\Subjects\Subject13\Timeseries\MRItimeseries\fMRItimeseries_S13_kw3.mat

Ignacio Rebollo 16/3/2015
commented 28/06/2017

%}

% kernelWidth = cfgMain.kernelWidth;

% Import

isrest = strcmp(cfgMain.task,'REST')
if isrest
dataDir = strcat(global_path2subject(subj_idx,cfgMain.sample),'fMRI',filesep); %navi
% dataDir = strcat(global_path2subject(subj_idx),'fMRI',filesep,'acquisition1',filesep,'RestingState',filesep); %phys

else
dataDir = strcat(global_path2subject(subj_idx,cfgMain.sample),'taskfMRI',filesep,cfgMain.task);
end

output_filename = global_filename(subj_idx,cfgMain,'BOLDTimeseriesFilename')

BOLDtimeseries = []; % "raw" fmri data
BOLDtimeseries.fsample  = 1/2;  %0.5 hz, ~ TR=2s


% identify the files corresponding to the deisred spatial kernel smoothing
% obtained from the output from mri preprocessing 


    filename= dir( fullfile( dataDir,'s3waf*.nii')); %# list all *.hdr files of preprocessed images, this are the 450 volumes
%     filename= dir( fullfile( dataDir,'s3waf*.img')); %# list all *.hdr files of preprocessed images, this are the 450 volumes


% list all files
filename = {filename(~[filename.isdir]).name}';

for i=1:length(filename)
  filename{i} = fullfile(dataDir, filename{i});% Only interested in the name of the files
end

% data has to be stored from a 3D matrix to one big vector tmp.anatomy(:,:,:)
datVector = zeros(length(filename),153594); % data in vector format(time,x*y*z) 
% 153594 correpond to the number of voxels of physiens EPI sequences (3mm)

for i=1:length(filename) %iterates through all volumes, load them, and put them in matrix
  disp(i);
  tmp = ft_read_mri(filename{i}); % calls fieldtrip to read the mri intensity images of each volume
  datVector(i,:) = tmp.anatomy (:); % concatenate them
end

BOLDtimeseries.time  = [0:2:(length(filename)*2)-1]; % create time axis
BOLDtimeseries.trialVector = datVector;
BOLDtimeseries.transform=tmp.transform(:,:);%transformation matrix; 

save(output_filename,'BOLDtimeseries')



% 
% 
% 
% % % list all files
% dataDir = strcat(global_path2subject(subj_idx),'rest3d',filesep);
% dataDirSubject= dir(dataDir) %# list all *.hdr files of preprocessed images, this are the 450 volumes
% filename = {dataDirSubject(~[dataDirSubject.isdir]).name}';
% 
% 
% 
% for i=1:length(filename)
%   filename{i} = fullfile(dataDir, filename{i});% Only interested in the name of the files
% end
% 
% % data has to be stored from a 3D matrix to one big vector tmp.anatomy(:,:,:)
% datVector = zeros(length(filename),153594); % data in vector format(time,x*y*z) 
% % 153594 correpond to the number of voxels of physiens EPI sequences (3mm)
% 
% for i=1:length(filename) %iterates through all volumes, load them, and put them in matrix
%   disp(i);
%   tmp = ft_read_mri(filename{i}); % calls fieldtrip to read the mri intensity images of each volume
%   datVector(i,:) = tmp.anatomy (:); % concatenate them
% end
% 
% BOLDtimeseries.time  = [0:2:length(BOLDVector)]; % create time axis
% BOLDtimeseries.trialVector = BOLDVector';
% BOLDtimeseries.transform=tmp.transform(:,:);%transformation matrix; 
% 
% save(output_filename,'BOLDtimeseries')
end