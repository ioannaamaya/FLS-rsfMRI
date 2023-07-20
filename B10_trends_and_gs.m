function B10_trends_and_gs(data_dir, filter, filtered_file, SPM_path, fbm)

f=spm_select('List',data_dir, filter);
rest_file=[data_dir filesep f];
[AllVolume, vsize, AllFileList ,Header, nVolumn] =rp_to4d(rest_file);
[nDim1, nDim2, nDim3, nDimTimePoints]=size(AllVolume);

%apply GM mask
[gm_mask, vsize, AllFileList ,Header, nVolumn] =rp_to4d(fbm);
data=[];
for i=1:nDimTimePoints
    temp=AllVolume(:,:,:,i);
    data(:,i)=temp(find(gm_mask)); % 2d, voxel x time
end


%check for scrubbing
n=1;
f3=spm_select('List',data_dir, ['^' filtered_file(n:end) '.*\_FWDstat.mat']);
while isempty(f3) && n<length(filtered_file)
    n=n+1;
    f3=spm_select('List',data_dir,['^' filtered_file(n:end) '.*\_FWDstat.mat']);
end
if ~isempty(f3)
    load([data_dir filesep f3],'outliers')
    %mask outliers
    data_masked=data(:,find(~outliers));
else
    data_masked=data;
end



%calculate global linear trend
fit_linear=polyfit(find(~outliers),squeeze(mean(data_masked))',1);
fit_quadratic=polyfit(find(~outliers),squeeze(mean(data_masked))',2);
trend_linear=polyval([fit_linear(1) fit_linear(2)],[1:nDimTimePoints])';
trend_quadratic=polyval([fit_quadratic(1) fit_quadratic(2) fit_quadratic(3)],[1:nDimTimePoints])';
%calculate global mean signal
gs=squeeze(mean(data));

%normalize regressors
trend_linear=(trend_linear-mean(trend_linear))/std(trend_linear);
trend_quadratic=(trend_quadratic-mean(trend_quadratic))/std(trend_quadratic);
gs=((gs-mean(gs))/std(gs))';

%save to disk
save([rest_file(1:end-4) '_trend_linear.txt'],'trend_linear','-ASCII','-DOUBLE','-TABS');
save([rest_file(1:end-4) '_trend_quadratic.txt'],'trend_quadratic','-ASCII','-DOUBLE','-TABS');
save([rest_file(1:end-4) '_global_signal.txt'],'gs','-ASCII','-DOUBLE','-TABS');




