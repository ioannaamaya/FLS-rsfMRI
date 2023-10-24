%%

clc
clear all
close all
addpath('E:\Data Analysis\ASFC_Data_processing')

src_dir      = 'E:\BIDS_Flicker\derivatives';
all_files = dir(fullfile(src_dir, 'sub*')); % looking for subject's folders
SJs = {all_files.name}; % creating cell array with subjects' IDs 


analysis_switch = [1:4];

condition_names_experimental = '*10Hz_run*bold.mat'; 

condition_names_placebo =  '*baseline_run*bold.mat'; 

exp_file_names2 = '10Hz_bold';
%names that will be given to the perc0 file. Use the name of experimental
%condition 
name = 'average_10Hz_bold_perc0';
contrast_name = 'R2Rconn150_10Hz_vs_baseline_'; 

%change here to use files with or without GSR 
prefix = 'R2Rconn150_Fh01l08_Rhclqg_s3'; 
prefix2 = '*Fh01l08_Rhclqg_s3*'; 
suffix = '*Fh01l08_Rhclqg_s3.mat'; 

    
%% Percent 0 
if ismember(1,analysis_switch)
    for subject = 1:numel(SJs)
        %display(['Step 13, cut data: ' cut_nifti_files(file).name])
        sub_dir = [all_files(subject).folder filesep all_files(subject).name];
        placebo_matrices =  dir(fullfile(sub_dir, '**', [prefix, condition_names_placebo])); 
        experimental_matrices = dir(fullfile(sub_dir, '**', [prefix, condition_names_experimental])); 
        currname = [SJs{subject} '_' name];
        E1_perc0_analysis(sub_dir, placebo_matrices, experimental_matrices, currname);
        clear currname
    end
end

%% Subtract matrices 
if ismember(2,analysis_switch)
    for subject = 1:numel(SJs)
        sub_dir = [all_files(subject).folder filesep all_files(subject).name];
        placebo_matrices =  dir(fullfile(sub_dir, '**', [prefix, condition_names_placebo])); 
        experimental_matrices = dir(fullfile(sub_dir, '**', [prefix, condition_names_experimental])); 
        name = [SJs{subject} '_' name];
        E2_subtract_matrices(sub_dir, placebo_matrices, experimental_matrices, exp_file_names2, SJs, subject, contrast_name);
    end
end

%% Average matrices 
if ismember(3,analysis_switch)
   all_matrices = dir(fullfile(src_dir, '**', [contrast_name, prefix2, 'mat'])); %search for the matrirces created in previous step
   E3_average_matrices(all_matrices, src_dir)
end

%% Visualzie matrix without old ROIs
if ismember(4,analysis_switch)
   E4_plot_conn_matrix(src_dir, contrast_name, suffix)
end
