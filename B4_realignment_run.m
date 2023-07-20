function B4_realignment_run(run_dir, filter)

warning off

fileset = {};
% select the files
f = spm_select('List', run_dir, filter);

% number of volumes
numVols = size(f,1);
% create SPM style file list for model specification
V=spm_vol([repmat([run_dir filesep], numVols, 1) f ',:']);

for i=1:size(V,1)
    fileset{i} = [repmat([run_dir filesep], numVols, 1) f repmat([',' int2str(i)] , numVols, 1)];
end
clear f;

matlabbatch{1}.spatial{1}.realign{1}.estwrite.data             = {fileset'};
matlabbatch{1}.spatial{1}.realign{1}.estwrite.eoptions.quality = 0.9;
matlabbatch{1}.spatial{1}.realign{1}.estwrite.eoptions.sep     = 4;
matlabbatch{1}.spatial{1}.realign{1}.estwrite.eoptions.fwhm    = 5;
matlabbatch{1}.spatial{1}.realign{1}.estwrite.eoptions.rtm     = 1;
matlabbatch{1}.spatial{1}.realign{1}.estwrite.eoptions.interp  = 2;
matlabbatch{1}.spatial{1}.realign{1}.estwrite.eoptions.wrap    = [0 0 0];
matlabbatch{1}.spatial{1}.realign{1}.estwrite.eoptions.weight  = {''};
matlabbatch{1}.spatial{1}.realign{1}.estwrite.roptions.which   = [2 1];
matlabbatch{1}.spatial{1}.realign{1}.estwrite.roptions.interp  = 4;
matlabbatch{1}.spatial{1}.realign{1}.estwrite.roptions.wrap    = [0 0 0];
matlabbatch{1}.spatial{1}.realign{1}.estwrite.roptions.mask    = 1;
matlabbatch{1}.spatial{1}.realign{1}.estwrite.roptions.prefix    = 'r';

inputs = cell(0, 1);

% save the job variable to disc
%fprintf(['Realignment ', SJ, '\n'])
spm('defaults', 'FMRI');
spm_jobman('serial', matlabbatch, '', inputs{:});
%spm_jobman('serial', matlabbatch, inputs{:});
clear jobs

