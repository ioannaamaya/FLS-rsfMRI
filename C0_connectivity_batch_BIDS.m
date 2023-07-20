% function C0_connectivity_batch_LUCIA_tn

% This is a Batch-Script to compute connectivity matrices for all subjects
% The pre-processed data for each run is taken individually
% connectivity matrices are saved in the current folder
% Alternatively one can use the CMs struct from the workspace for plotting
% and further processing
%Toolbox required: Bramilla

%TTS/TN: wurde irgendwie von superbatch ?bergeben
%atlas = fullfile(data_base, 'V1_mask_Carhart_Harris.nii');
%data_base ='1st_level_V1';
%regressor_prefix ='V1_ROI_time_course_'

clc
clear all
close all
addpath('E:\Data Analysis\ASFC_Data_processing')
addpath(genpath('E:\Toolboxes\bramila-master'))
%#####################################################
%#################### INPUT ##########################
%#####################################################

%SPM-path
SPM_path  = addpath('E:\Toolboxes\spm12');

%data source directory
src_dir      = 'E:\BIDS_Flicker\derivatives';
%segment = '_end_';
nifti_files = dir(fullfile(src_dir, '**', ['Fh01l08_Rhclq_s3', '*bold.nii']));%look for all processed nifti files


% specify full path to the ROI NIFTI file (atlas)
atlas = 'E:\Data Analysis\Atlases\atlas_finale_3.nii';
%define grey matter mask
mask= 'E:\Data Analysis\brain_mask.nii';

% selection of analysis steps to be performed
analysis_switch = [2];

%# step 1  reslice atlas with ROIs
%  choose reference subject and run (used for reslicing)
ref_sub=1;
ref_run=1;

    
%# step 2  ROI2ROI, calculate between-ROI functional connectivity
%  specify which ROIs of the atlas to include
ROI_values = [1:150]';

%# step 3  seed-based analysis (for full-brain seed-based connectivity)
%  specify seed-region based on atlas-regions or file
%
%  e.g. for labels_Neuromorphometrics.nii
%  parahippocampus:  seed_ROI=[170, 171];
%  calcarine sulcus / V1: seed_ROI=[108, 109];
%  vmPFC: seed_ROI=[146,147,178,179,104,105,124,125,136,137];
%  Thalamus r/l: seed_ROI=[92 93] 
% seed_ROI=[93];
%
%  or specify file: e.g. seed_ROI='H:\Ganzfeld\Data\ECM_flexFact_fastECM_Fh01l08_Rhclq_sm0.4wrax3f4d\seed_Pcun_ECM_pre_vs_ganz.nii';
% seed_ROI='H:\Ganzfeld\Data\ECM_flexFact_fastECM_Fh01l08_Rhclq_sm0.4wrax3f4d\seed_conj_lePPC.nii';

%# step 4  calculate ECMs (using fastECM)
ztransform=0;   % 1=yes / 0=no
smooth=1;       % 1=yes / 0=no

%#####################################################
%#################### INPUT end ######################
%#####################################################

%% reslice the atlas to match functional data, if necessary
if ismember(1,analysis_switch)
    sub_run_path = [nifti_files(ref_sub).folder];
    f = spm_select('List', sub_run_path, ['^' nifti_files(ref_sub).name]);
    rest_path = [sub_run_path filesep f];
    C1_check_reslice(rest_path, atlas);
end
%%
%  atlas = 'E:\Data Analysis\Atlases\resliced_AAL3v1.nii';
atlas = 'E:\Data Analysis\Atlases\hybrid_atlas.nii'
%% ROI2ROI functional connectivity
if ismember(2,analysis_switch)
    C2_ROI2ROI_conn_masked(nifti_files, src_dir,atlas,ROI_values);
end

%% seed-based analysis
if ismember(3,analysis_switch)
    C3_ROI2voxel_conn_masked(nifti_files,src_dir,atlas,seed_ROI,mask)
end

%% calculate ECM maps
if ismember(4,analysis_switch)
    C4_fast_ecm(nifti_files,mask,ztransform,smooth)
end

