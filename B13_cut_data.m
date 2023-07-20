%Cut data in 5min segments at the beginning, end, and middle of the
%scanning 
function B13_cut_data(src_dir, TR, nifti_file, segment_length, segment_start)

segment_size = round(segment_length/TR); %5 minutes segment 

%nifti_files = dir(fullfile(src_dir, '**', ['Fh01l08_Rhclq_s6m0.4wrasub-', '*bold.nii']));

    tempdir= src_dir;
    cd(tempdir)
    f = spm_select('List',tempdir, ['^' nifti_file]);
    nr_scans=length(spm_vol([tempdir filesep f]));
    current_TR=segment_start;
    %cut first segment of the data
    if (current_TR + segment_size-1) <= nr_scans
        display(['file ' nifti_file ', segment ' num2str(current_TR) '-' num2str(current_TR + segment_size-1)])
        M = load_untouch_nii([tempdir filesep f],[current_TR:current_TR + segment_size-1]);
        M.fileprefix = [M.fileprefix '_start_TR_' num2str(current_TR) '_' num2str(current_TR + segment_size-1)];
        save_untouch_nii(M,[M.fileprefix '.nii']);
    end
    %cut last segment of the data
    current_TR = (nr_scans - segment_size) + 1; 
    if current_TR ~= segment_start %if data is <5 min, current_TR should be 1 and cutting should be skipped
        display(['file ' nifti_file ', segment ' num2str(current_TR) '-' num2str(current_TR + segment_size-1)])
        M = load_untouch_nii([tempdir filesep f],[current_TR:nr_scans]);
        M.fileprefix = [M.fileprefix '_end_TR_' num2str(current_TR) '_' num2str(current_TR + segment_size-1)];
        save_untouch_nii(M,[M.fileprefix '.nii']);
    end
    %cut middle 5 min of the data
    current_TR = (round(nr_scans/2)) - round(segment_size/2) + 1;
    if current_TR ~= segment_start %if data is <5 min, current_TR should be 1 cutting should be skipped
        display(['file ' nifti_file ', segment ' num2str(current_TR) '-' num2str(current_TR + segment_size-1)])
        M = load_untouch_nii([tempdir filesep f],[current_TR:current_TR + segment_size-1]);
        M.fileprefix = [M.fileprefix '_middle_TR_' num2str(current_TR) '_' num2str(current_TR + segment_size-1)];
        save_untouch_nii(M,[M.fileprefix '.nii']);
    end
    

