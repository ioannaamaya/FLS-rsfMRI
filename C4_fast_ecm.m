function C4_fast_ecm(nifti_files,mask,ztransform,smooth)

% This is a Batch-Script to compute Eigen-Centrality-Maps for all subjects
% The pre-processed data for each run is taken individually

mask_grey=load_untouch_nii(mask);
gg=mask_grey.img;

for file=1:length(nifti_files)
    sub_dir=[nifti_files(file).folder];
    ECM_dir=[sub_dir filesep 'ECM_results'];
    mkdir(ECM_dir)
    filter =  (extractBefore(nifti_files(file).name, "_start"));    
    filter2 = [nifti_files(file).name];
    filter3 =  (extractBefore(nifti_files(file).name, ".nii"));
    scan_indx=[];
    %check for scrubbing
    n=1;
    f2=spm_select('List',sub_dir,['^' filter(n:end) '.*\_FWDstat.mat']);
    while isempty(f2) && n<length(filter)
        n=n+1;
        f2=spm_select('List',sub_dir,['^' filter(n:end) '.*\_FWDstat.mat']);
    end
    if ~isempty(f2)
        load([sub_dir filesep f2],'outliers')
    end
    %find functional files
    f1 = spm_select('List',sub_dir, ['^' filter2]);
    for f=1:size(f1,1)
        tempfile=deblank(f1(f,:));
        display([file ', file ' tempfile]);
        functional_SJdata=[sub_dir filesep tempfile];
            
        %check data segment
        n=1;
        while ~strcmp('TR',tempfile(n:n+1)) && n<length(tempfile)-1
            n=n+1;
        end
        if n==length(tempfile)-1
            inx=[1 length(spm_vol(functional_SJdata))];
        else
            while ~strcmp('_',tempfile(n))
                n=n+1;
            end
                
            tmp=[];
            n=n+1;
            while ~strcmp('_',tempfile(n))
                tmp=[tmp tempfile(n)];
                n=n+1;
            end
            inx(1)=str2num(tmp);
                
            tmp=[];
            n=n+1;
            while ~strcmp('.',tempfile(n))
                tmp=[tmp tempfile(n)];
                n=n+1;
            end
            inx(2)=str2num(tmp);
        end
        if ~isempty(f2)
            scan_indx=find(~outliers(inx(1):inx(2)));
        end
            
        fastECM_hacked(functional_SJdata,0,1,0,24,mask,0,scan_indx);
            
         %### Find ECM files for renaming
         files2move = cellstr(spm_select('List',sub_dir, [ '^*' filter3 '.*\ECM.nii'] ));
         for i=1:length(files2move)
             movefile(fullfile(sub_dir, files2move{i}),[sub_dir filesep files2move{i}(end-10:end-4) '_' files2move{i}(1:end-12) '.nii']);
         end
            
         %###  Z-Transform and smoothing
         file2transform=spm_select('List',sub_dir,['^fastECM_.*\' filter3]);
            
         if ztransform
             temp=load_untouch_nii([sub_dir filesep file2transform]);
             tt=temp.img(find(gg));
             global_mean=mean(tt);
             global_std=std(tt);
             temp.img=(temp.img - global_mean) / global_std;
             temp.img(find(gg==0))=0;
             temp.fileprefix=[sub_dir filesep 'z_' file2transform];
             save_untouch_nii(temp,temp.fileprefix);
             file2transform=['z_' file2transform];
         end
            
         if smooth
             %###  Smoothing
             matlabbatch{1}.spm.spatial.smooth.data = {[sub_dir filesep file2transform]};
             matlabbatch{1}.spm.spatial.smooth.fwhm = [6 6 6];
             matlabbatch{1}.spm.spatial.smooth.dtype = 0;
             matlabbatch{1}.spm.spatial.smooth.im = 0;
             matlabbatch{1}.spm.spatial.smooth.prefix = 's';
                
             spm('defaults', 'FMRI');
             spm_jobman('serial', matlabbatch);
             clear matlabbatch
         end
         movefile([sub_dir filesep '*ECM*.nii'], [ECM_dir]) 
    end
end
