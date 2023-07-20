function prefix = B11_regress_out_nuisance(data_dir, filter, fbm, hm, cc, tl, tq, gs)

% make NIFTIs with regressed out nuisance variables

%% load covariates
%head motion
if hm
    dinfo='> head motion ';
    prefix='Rh';
    f=spm_select('List',data_dir,['^rp_' filter '.*\.txt']);
    n=1;
    while isempty(f)
        f=spm_select('List',data_dir,['^rp_' filter(n:end) '.*\.txt']);
        n=n+1;
    end
    if ~isempty(f)
        nCOV=load([data_dir filesep f]);
    else
        nCOV=[];
    end
else
    nCOV=[];    %nuisance Covariates
    dinfo='';
    prefix='R';
end

%CompCorr
if cc
    f=spm_select('List',data_dir,['^' filter '.*\_CompCorPCs.txt']);
    n=1;
    while isempty(f)
        f=spm_select('List',data_dir,['^' filter(n:end) '.*\_CompCorPCs.txt']);
        n=n+1;
    end
    if ~isempty(f)
        temp = load([data_dir filesep f]);
        nCOV=[nCOV temp];
    end
    clear temp
    dinfo=[dinfo ' > CompCorr PCs '];
    prefix=[prefix 'c'];
end

%linear trend
if tl
    f=spm_select('List',data_dir,['^' filter '.*\_trend_linear.txt']);
    temp = load([data_dir filesep f]);
    nCOV=[nCOV temp];
    clear temp
    dinfo=[dinfo ' > linear trend '];
    prefix=[prefix 'l'];
end

%quadratic trend
if tq
    f=spm_select('List',data_dir,['^' filter '.*\_trend_quadratic.txt']);
    temp = load([data_dir filesep f]);
    nCOV=[nCOV temp];
    clear temp
    dinfo=[dinfo ' > quadratic trend '];
    prefix=[prefix 'q'];
end

%global signal
if gs
    f=spm_select('List',data_dir,['^' filter '.*\_global_signal.txt']);
    temp = load([data_dir filesep f]);
    nCOV=[nCOV temp];
    clear temp
    dinfo=[dinfo ' > global signal '];
    prefix=[prefix 'g'];
end

prefix=[prefix '_'];
%% load data

f=spm_select('List',data_dir,['^' filter '.*\.nii']);
rest_file=[data_dir filesep f];
[AllVolume, vsize, AllFileList ,Header, numVols] =rp_to4d(rest_file);

%load full brain mask
[gm_mask, vsize, AllFileList ,Header, nVolumn] =rp_to4d(fbm);
[nDim1, nDim2, nDim3]=size(gm_mask);

%% check for scrubbing
n=1;
f3=spm_select('List',data_dir,['^' filter(n:end) '.*\_FWDstat.mat']);
while isempty(f3) && n<length(filter)
    n=n+1;
    f3=spm_select('List',data_dir,['^' filter(n:end) '.*\_FWDstat.mat']);
end
if ~isempty(f3)
    load([data_dir filesep f3])
    %mask 'bad' time points (outliers)
    nCOV_masked = nCOV(find(~outliers),:);
    AllVolume_masked = AllVolume(:,:,:,find(~outliers));
    numVols_masked=size(AllVolume_masked,4);
else
    nCOV_masked = nCOV;
    AllVolume_masked = AllVolume;
    numVols_masked=size(AllVolume_masked,4);
end


%% regress out nuisance covariates
fprintf('\n');
display(['regress out nuisance covariates: ' dinfo]);

AllVolume=reshape(AllVolume,[],numVols)';    % Convert into 2D
AllVolume_masked=reshape(AllVolume_masked,[],numVols_masked)';    % Convert into 2D


Residual_data=[];
for i = 1:size(AllVolume_masked,2)
    resp = AllVolume_masked(:,i);
    s = regstats(resp,nCOV_masked,'linear',{'beta'});
    %Betas = s.beta(1:[length(s.beta)-1]);
    Betas = s.beta(2:(size(nCOV,2)+1));
    nuisance = nCOV*Betas;
    Residual_data(:,i) = AllVolume(:,i) - nuisance;
end

%mask out non-brain regions
Residual_data(find(~gm_mask))=0;
% Convert into 4D
AllVolumeBrain=reshape(Residual_data',[nDim1, nDim2, nDim3, numVols]);
% Save all images to disk
fprintf('\n\t Saving data.\tWait...');
Header.pinfo = [1;0;0];
rp_Write4DNIfTI(AllVolumeBrain,Header,[data_dir filesep prefix f]);
%clear All* temp_path

fprintf('...done\n\n');