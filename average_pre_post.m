%% average pre & post 

clc
clear all
close all
src_dir = 'E:\BIDS_Flicker\derivatives';

all_files = dir(fullfile(src_dir, 'sub*')); % looking for subject's folders
SJs = {all_files.name}; % creating cell array with subjects' IDs 

% subtract matrices
for subject = 1:numel(SJs)
    sub_dir = [all_files(subject).folder filesep all_files(subject).name];
    pre_matrix =  dir(fullfile(sub_dir, '**', ['R2Rconn150_Fh01l08_Rhclqg_s3', '*pre*bold.mat'])); 
    post_matrix =  dir(fullfile(sub_dir, '**', ['R2Rconn150_Fh01l08_Rhclqg_s3', '*post*bold.mat']));
    A = [pre_matrix.folder filesep pre_matrix.name];
    matrix_A = load(A);
    B = [post_matrix.folder filesep post_matrix.name];
    matrix_B = load(B);
    CorrMat =  ((matrix_B.CorrMat + matrix_A.CorrMat)/2);
    ROI = matrix_A.ROI;
    Percent0 = matrix_A.percent0;
    new_name = (extractAfter(B, 'hybrid_atlas\'));
    new_name = strrep(new_name, 'task-post' , 'task-baseline');
    eval(['save ' pre_matrix.folder filesep new_name ' CorrMat ROI Percent0']);
end 
