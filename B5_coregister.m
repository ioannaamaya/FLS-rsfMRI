function B5_coregister(func_dir, struct_dir, filter_mean, filter_struct, func_files)

warning off

f1 = spm_select('List', struct_dir, filter_struct);
numVols = size(f1,1);
structural = cellstr([repmat([struct_dir filesep], numVols, 1) f1 repmat(',1', numVols, 1)]);

f2 = spm_select('List', func_dir, filter_mean);
numVols = size(f2,1);
mean_img   = cellstr([repmat([func_dir filesep], numVols, 1) f2 repmat(',1', numVols, 1)]);

Images=cellstr(spm_select('ExtFPList', func_dir,func_files,Inf));

%--------------------------------------------------------------------------
% --------------------------- Coregister       ----------------------------
%--------------------------------------------------------------------------
matlabbatch{1}.spm.spatial.coreg.estimate.ref = structural;
matlabbatch{1}.spm.spatial.coreg.estimate.source = mean_img;
matlabbatch{1}.spm.spatial.coreg.estimate.other = Images;
matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.cost_fun = 'nmi';
matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.sep = [4 2];
matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.tol = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.fwhm = [7 7];

% run job
spm_jobman('run', matlabbatch)
clear jobs