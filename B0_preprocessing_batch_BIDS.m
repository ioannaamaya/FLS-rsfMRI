% function B0_preprocessing_batch_BIDS

% please send questions to Till Nierhaus (till.nierhaus@fu-berlin.de) or Timo T. Schmidt ().

%### step A, structure data before running the preprocessing
%### --> convert DICOM images to 4D NIFTI image
%### --> subject folders (e.g. 'sub-01') including sessions (e.g. 'ses-1'), 
%functional (e.g. 'func') and anatomy (e.g. 'anat')
%############################################################################################

% required toolboxes:
% SPM12 ( http://www.fil.ion.ucl.ac.uk/spm/ )
% Rest_plus ( http://restfmri.net/forum/index.php?q=rest )
% hMRI (https://www.cbs.mpg.de/departments/neurophysics/software/hmri-toolbox)
% Step 7  Scrubbing: BRAMILA toolbox (  )
clc
clear all
close all
addpath('E:\Data Analysis\ASFC_Data_processing')
addpath(genpath('E:\Toolboxes\hMRI-toolbox-master'))
addpath(genpath('E:\Toolboxes\bramila-master'))
addpath(genpath('E:\Toolboxes\210630_1608_RESTplus_v1.25'))
addpath('E:\Toolboxes\spm12')

%#####################################################
%#################### INPUT ##########################
%#####################################################

%SPM-path
SPM_path  ='E:\Toolboxes\spm12';

%data source directory
src_dir      = 'E:\BIDS_Flicker\derivatives';        

all_files = dir(fullfile(src_dir, 'sub*')); % looking for subject's folders
SJs = {all_files.name}; % creating cell array with subjects' IDs 

nifti_files = dir(fullfile(src_dir, '**', ['sub-', '*bold.nii'])); %look for all functional nifti files
anat_files = dir(fullfile(src_dir, '**', ['sub-', '*T1w.nii'])); %look for all anat nifti files

%anatomical masks
wm_mask=['E:\Data Analysis\wm_mask_eroded.nii']; %white matter mask file
csf_mask=['E:\Data Analysis\csf_mask_eroded.nii']; %csf mask file
full_brain_mask=['E:\Data Analysis\full_brain_mask.nii']; %full brain mask file

% selection of analysis steps (1-12) to be performed
analysis_switch = [1:12];
start_prefix = 'sub'; 

%now we get the data from the json file 
json_files = (dir(fullfile(src_dir, '**', ['task', '*json']))); %extract all json files, althoguh they should have the same info 
%because for some datasets (flicker and Ganzfeld) we have json files for
%each nift file and these are named differently, we have to check if the
%first command returns an empty structure. If yes, it means the json files
%have a different naming, starting with subject 
if isequal(size(json_files), [0, 1]) 
    json_files = (dir(fullfile(src_dir, '**', ['sub-', '*bold.json'])));
end

json_file = [json_files(1).folder, filesep, json_files(1).name]; %we select the first json file to extract metadata from 
TR_json = get_metadata_val(json_file,'RepetitionTime') / 1000; % repetition time in sec
slice_timing = get_metadata_val(json_file,'SliceTiming'); %extract slice timing
n_slices_json = height(slice_timing); %compute number of slices from slice timing
[~,y]= sort(slice_timing); %compute slice order 
slice_order = y';

%now get the same info from nifti header  
nifti_file_metadata = [nifti_files(1).folder, filesep, nifti_files(1).name]; 
info = niftiinfo(nifti_file_metadata);
TR_nifti = info.PixelDimensions(4); 
n_slices_nifti = info.ImageSize(3);

%compare json and nifti header 
if round(TR_nifti, 4) ~= round(TR_json, 4)
    warning ("TR does not match between json file and nifti") 
end 
if n_slices_json ~= n_slices_nifti
    warning ("Number of slices does not match between json file and nifti") 
end 

%%% I suggest using the json values at least for TR since we know it is
%%% missing in the nifti headers for some datasets
TR = TR_json;
n_slices = n_slices_json; 

%# step 1  Segmentation
%# ------ Create nuisance masks on your own or take the provided ones
% create gray and white matter images and bias-field corrected structural image
%# ------ NOTE: spike removal (e.g. "artrepair") should be performed as first step
%# step 2 --> remove first x scans                       --> prefix: x(number of cut volumes)
    x=0;
%# step 3 --> slice time correction                      --> prefix: a
%  we first do slice time correction, then realignment, because of interleaved slice order
%Correct differences in image acquisition time between slices 
    refslice=slice_order(round(length(slice_order)/2)); % reference slice, middle slice in slice_order
%# step 4  Realignment                                --> prefix: r
% estimate the 6 parameter (rigid body) spatial transformation
%that will align the times series of images and will modify the header of the input images
%(*.hdr), such that they reflect the relative orientation of the data after correction for movement
%artefacts.
%# step 5  Coregister (estimate) mean-epi 2 anatomy
% implement a coregistration between the structural and functional data that 
%maximises the mutual information
%# step 6  Normalization                              --> prefix: w
    vox_size=[2 2 2]; % voxel size in mm
%# step 7  Scrubbing: calculate, interpolate outliers --> prefix: m(scrub_thresh)
    scrub_thresh=0.4; % threshhold FD for scrubbing
%# step 8 Calculate WM and CSF Nuisance Signal
    numComp = 5; % number of principle components    
%# step 9 Smoothing                                  --> prefix: s
    kernel_size=[3 3 3]; %FWHM kernel size
%# step 10 Calculate Trends & Global Signal
%# step 11 Regress out nuisance covariates
    hm=1;   % head motion parameters from realignment (step 4)
    cc=1;   % CompCorr WM and CSF principal components (step 8)
    tl=1;   % linear trend (step 10)
    tq=1;   % quadratic trend (step 10)
    gs=0;   % global signal (step 10)
%# step 12 Bandpass Filter                            --> prefix: hp()_lp()
    hpf = 0.01; %LowCutoff / high pass filter in Hz
    lpf = 0.08; %HighCutoff / low pass filter in Hz
%# Step 13 Cut data: cut three segmets, first/middle/last
    segment_length = 300; %5 min segment in seconds
    segment_start = 1; %index of the first TR where the first segment starts

%#####################################################
%#################### INPUT end ######################
%#####################################################

%%
currPrefix=start_prefix;
currPrefix = ['s3m0.4wrasub']
%% Segmentation
if ismember(1,analysis_switch)
    warning off
    for file = 1:numel(anat_files) 
        display(['Step 1, segmentation: ' anat_files(file).name])
        struct_dir = [anat_files(file).folder];
        B1_segmentation(struct_dir, SPM_path, ['^' anat_files(file).name]);
    end
end

%% Delete first X scans
if ismember(2,analysis_switch)
    for file = 1:numel(nifti_files)
        display(['Step 2, delete first ' num2str(x) ' volumes: ' nifti_files(file).name])
        run_dir = [nifti_files(file).folder];
        B2_delete_scans(run_dir, ['^' nifti_files(file).name], 2);
    end
    currPrefix=['x' num2str(x) currPrefix];
end

%% Slice time correction
a_nifti_files = dir(fullfile(src_dir, '**', [currPrefix, '*bold.nii']));
if ismember(3,analysis_switch)
    for file = 1:numel(a_nifti_files) 
        display(['Step 3, slice time correction: ' a_nifti_files(file).name])
        %file.task = strcat(SJs{sj}, '_task-', tasks{r}, '_bold.nii');
        run_dir = [a_nifti_files(file).folder];
        B3_slice_time_correction(run_dir, ['^' a_nifti_files(file).name], n_slices,slice_order,refslice,TR);
    end
    currPrefix=['a' currPrefix];
end

%% Realignment
r_nifti_files = dir(fullfile(src_dir, '**', [currPrefix, '*bold.nii']));
if ismember(4,analysis_switch)
    for file = 1:numel(r_nifti_files) 
        display(['Step 4, realignment: ' r_nifti_files(file).name])
        run_dir = [r_nifti_files(file).folder];
        B4_realignment_run(run_dir, ['^' r_nifti_files(file).name]);
    end
    meanPreffix = ['mean' currPrefix];
    rpPrefix = ['rp_' currPrefix];
    currPrefix=['r' currPrefix];
end

%% Coregister (estimate) mean-epi 2 anatomy
if ismember(5,analysis_switch)
    meanPreffix = 'meanasub';
    warning off
    for anat = 1:numel(anat_files) 
        c_nifti_files = dir(fullfile(src_dir, '**', [currPrefix SJs{anat}(4:end) '*bold.nii']));
        mean_nifti_file = dir(fullfile(src_dir, '**', [meanPreffix SJs{anat}(4:end) '*bold.nii']));
        for file = 1:numel(c_nifti_files) 
            func_dir = [c_nifti_files(file).folder];
            struct_dir = [anat_files(anat).folder];
            B5_coregister(func_dir, struct_dir, ['^' mean_nifti_file(file).name], ['^' anat_files(anat).name], ['^' c_nifti_files(file).name]);
        end
    end
end

%% Normalization
if ismember(6,analysis_switch)
    for anat = 1:numel(anat_files) 
        w_nifti_files = dir(fullfile(src_dir, '**', [currPrefix SJs{anat}(4:end) '*bold.nii']));
        struct_dir = [anat_files(anat).folder];
        for file = 1:numel(w_nifti_files) 
            display(['Step 6, normalization: ' w_nifti_files(file).name])
            data_dir = [w_nifti_files(file).folder];
            B6_normalization_run(data_dir, struct_dir, ['^' w_nifti_files(file).name], vox_size);
        end
    end
    currPrefix=['w' currPrefix];
end

%% Scrubbing: calculate outliers
m_nifti_files = dir(fullfile(src_dir, '**', [currPrefix, '*bold.nii']));
rpPrefix = 'rp_asub';
if ismember(7,analysis_switch)
    rp_files = dir(fullfile(src_dir, '**', [rpPrefix, '*bold.txt']));
    scrub_prefix=['m' num2str(scrub_thresh)];
    for file = 1:numel(m_nifti_files)
        display(['Step 7, scrubbing: ' m_nifti_files(file).name])
        data_dir = [m_nifti_files(file).folder];
        %estimate and save motion statistics
        n=1;
        %%rp_file = (extractBefore(rp_files(file).name, ".nii"))
        f=spm_select('List', data_dir, ['^' rp_files(file).name]);
        while isempty(f)
            n=n+1;
            f=spm_select('List', data_dir, ['^' rp_files(file).name]);
        end
        cfg.motionparam=[data_dir filesep f];
        cfg.prepro_suite = 'spm';
        [fwd,rms]=bramila_framewiseDisplacement(cfg);
        outliers=fwd>scrub_thresh;
        percent_out=(sum(outliers)/length(outliers))*100;
        disp(['outliers for ' num2str(m_nifti_files(file).name) ': ' num2str(percent_out) '%']);
        filtered_file = (extractBefore(m_nifti_files(file).name, ".nii"));
        save([data_dir filesep scrub_prefix filtered_file '_FWDstat.mat'],'fwd','rms','outliers','percent_out','scrub_thresh','cfg')
        %srub outliers by replacing them with average of nearest neighbors
        B7_scrub_data(data_dir, ['^' m_nifti_files(file).name], outliers,  scrub_prefix);
        all_percent_out(file)=percent_out;
        all_rp{file}=load(cfg.motionparam);
        
    end
    currPrefix=[scrub_prefix currPrefix];
    save([src_dir filesep 'all_MOTIONstat_' currPrefix '.mat'],'m_nifti_files','scrub_thresh','all_percent_out','all_rp')
end

%% CompCorr
comp_nifti_files = dir(fullfile(src_dir, '**', [currPrefix, '*bold.nii']));
if ismember(8,analysis_switch)
    for file = 1:numel(comp_nifti_files)
        display(['Step 8, CompCorr: ' comp_nifti_files(file).name])
        data_dir = [comp_nifti_files(file).folder];
        B8_compcorr_run(data_dir, ['^' comp_nifti_files(file).name], numComp, wm_mask, csf_mask);
    end
end

%% Smoothing
s_nifti_files = dir(fullfile(src_dir, '**', [currPrefix, '*bold.nii']));
if ismember(9,analysis_switch)
    for file = 1:numel(s_nifti_files)
        display(['Step 9, smoothing: ' s_nifti_files(file).name])
        run_dir = [s_nifti_files(file).folder];
        B9_smoothing_run(run_dir, ['^' s_nifti_files(file).name], kernel_size);
        display([s_nifti_files(file).name ' is done'])
    end
    currPrefix=['s' num2str(unique(kernel_size)) currPrefix];
end
%% Calculate trends and global signal
trends_nifti_files = dir(fullfile(src_dir, '**', [currPrefix, '*bold.nii']));
if ismember(10,analysis_switch)
    for file = 1:numel(trends_nifti_files)
        display(['Step 10, calculate trends and global mean: ' trends_nifti_files(file).name])
        data_dir = [trends_nifti_files(file).folder];
        filtered_file = (extractBefore(trends_nifti_files(file).name, ".nii"));
        B10_trends_and_gs(data_dir, ['^' trends_nifti_files(file).name], ['^' filtered_file], SPM_path, full_brain_mask);
        %display([trends_nifti_files(file).name ' is done'])
    end
end

%% Regress out nuisance covariates
reg_nifti_files = dir(fullfile(src_dir, '**', [currPrefix, '*bold.nii']));
if ismember(11,analysis_switch)
    for file = 1:numel(reg_nifti_files)
        display(['Step 11, regress out nusiance covariates: ' reg_nifti_files(file).name])
        data_dir = [reg_nifti_files(file).folder];
        filtered_file = (extractBefore(reg_nifti_files(file).name, ".nii"));
        prefix = B11_regress_out_nuisance(data_dir, filtered_file, full_brain_mask, hm, cc, tl, tq, gs);
        
    end
    currPrefix=[prefix currPrefix];
end

%% Filter
filter_nifti_files = dir(fullfile(src_dir, '**', [currPrefix, '*bold.nii']));
if ismember(12,analysis_switch)
    for file = 1:numel(filter_nifti_files)
        display(['Step 12, filter data: ' filter_nifti_files(file).name])
        data_dir = [filter_nifti_files(file).folder];
        B12_bandpass_filter_run(data_dir, TR, hpf, lpf, ['^' filter_nifti_files(file).name]);
    end
    currPrefix = ['Fh01l08_' currPrefix];
end

%% Cut data
cut_nifti_files = dir(fullfile(src_dir, '**', [currPrefix, '*bold.nii']));
if ismember(13,analysis_switch)
    for file = 1:numel(cut_nifti_files)
        display(['Step 13, cut data: ' cut_nifti_files(file).name])
        data_dir = [filter_nifti_files(file).folder];
        B13_cut_data(data_dir, TR, ['^' cut_nifti_files(file).name], segment_length, segment_start);
    end
end
%%
% adjust the following to delete files with certain prefixes from the subject
% folders
if ismember(99,analysis_switch)
    for file = 1:numel(nifti_files)
        data_dir = [nifti_files(file).folder];
        cd(data_dir)
        %delete('Rhclq_s6m0.4wrasub*start_TR_1_100.nii')
        %delete('*start_TR_1_100.nii')
        %delete('rp_sub*.txt')
        %delete('d*.mat')
        %delete('s*.mat')
        %delete('w*.mat')     
    end
end