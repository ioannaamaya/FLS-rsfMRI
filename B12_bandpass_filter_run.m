function B12_bandpass_filter_run(data_dir, TR, hpf, lpf, filter_imgs)

spm('defaults','fmri');
spm_jobman('initcfg');

warning off
f = spm_select('List', data_dir, filter_imgs);
numVols = size(f,1);
Images=cellstr([repmat([data_dir filesep], numVols, 1) f]);

%##### 5) Temporal Filter
%display([temps ': Temporal Filter'])
%fprintf(['\n\t ' SJ ' Band Pass Filter working.\tWait...']);

%[AllVolume, vsize, AllFileList ,Header, nVolumn] =rp_to4d([data_dir filesep 'regressed.nii']);
[AllVolume, vsize, AllFileList ,Header, nVolumn] =rp_to4d(Images);
[nDim1, nDim2, nDim3, nDimTimePoints]=size(AllVolume);
% Convert into 2D
AllVolume=reshape(AllVolume,[],nDimTimePoints)';
%Remove the mean
AllMean=mean(AllVolume);
AllVolume=AllVolume-repmat(AllMean,[nDimTimePoints,1]);
%Filter
AllVolume = rp_IdealFilter(AllVolume, TR, [hpf, lpf]);  % describe in thesis %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Add the mean back after filter.
AllVolume=AllVolume+repmat(AllMean,[nDimTimePoints,1]);
% Convert into 4D
%AllVolumeBrain=reshape(int16(AllVolume)',[nDim1, nDim2, nDim3, nDimTimePoints]);
AllVolumeBrain=reshape(AllVolume',[nDim1, nDim2, nDim3, nDimTimePoints]);
%Save all images to disk
%fprintf(['\n\t ' SJ ' Saving filtered images.\tWait...']);
%Header.pinfo = [1;0;0];

high_str=num2str(lpf);
low_str=num2str(hpf);
prefix=['F' 'h' low_str(3:end) 'l' high_str(3:end) '_'];

rp_Write4DNIfTI(AllVolumeBrain,Header,[data_dir filesep prefix f]);
display('done');
clear Images Head* All* n* vsize
end