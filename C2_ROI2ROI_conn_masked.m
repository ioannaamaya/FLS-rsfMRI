function C2_ROI2ROI_conn_masked(functional_files, src_dir,atlas,ROI_values)

%%%%%%%%%%%%%
% with temporal masking of outliers (scrubbing)

% load atlas
mask_atl=load_untouch_nii(atlas);
atl=mask_atl.img;
[atlas_dir atlas_file atlas_ext]=fileparts(atlas);

%specify ROIs
field='ID';
ROI = struct(field, ROI_values);

% LOOP over files
for file=1:length(functional_files)
    sub_dir=[functional_files(file).folder];
    filter =  (extractBefore(functional_files(file).name, "_bold"));
    filter = [filter '_bold'];
    filter2 = [functional_files(file).name];
    outputdir=[sub_dir filesep 'ROI2ROI_FC_' atlas_file]; 
    mkdir(outputdir);
        
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
      % list functional data files
    f1 = spm_select('List', sub_dir, ['^' filter2]);
    for f=1:size(f1,1)
        tempfile=deblank(f1(f,:));
        display([file ', file ' tempfile]);
        % Loading preprocessed functional DATA
        fwdata=load_untouch_nii([sub_dir filesep tempfile]);
        %check data segment
        n=1;
        while ~strcmp('TR',tempfile(n:n+1)) && n<length(tempfile)-1
            n=n+1;
        end
        if n==length(tempfile)-1
            inx=[1 size(fwdata.img,4)];
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
            display(['scrub data ' num2str(inx(1)) ' - ' num2str(inx(2))])
           % apply temporal mask to the data
           fwdata.img(:,:,:,find(outliers(inx(1):inx(2))))=[];
        else
         display('#####################################################')
         display('########## No FWD file found for scrubbing ##########')
         display('#####################################################')
        end
        [nDim1, nDim2, nDim3, nDimTimePoints]=size(fwdata.img);
         % Convert into 2D
         AllVolume=reshape(fwdata.img,[],nDimTimePoints)';
         %extract ROI-mean signals
         AvgMat=zeros(length(ROI.ID),nDimTimePoints);
         for i=1:length(ROI.ID)
             temp_idx=find(atl==ROI(1).ID(i));
             x1=AllVolume(1,temp_idx);   
             percent0(i)=length(find(round(x1)==0))/length(x1);
             AvgMat(i,:) = mean(AllVolume(:,temp_idx),2);
         end
         %correlate ROI-mean signals
         CorrMat=zeros(length(ROI.ID),length(ROI.ID));
         for i=1:length(ROI.ID)
             temp=repmat(AvgMat(i,:),length(ROI.ID),1);
             CorrMat(:,i)=x_calcPairCorr(AvgMat,temp);
         end

         % Round (to avoid output of complex numbers) and apply Fisher z-transformation to all ROI-to-ROI pairs.
         CorrMat = atanh(round(CorrMat, 7));
            
         eval(['save ' outputdir filesep 'R2Rconn150_' tempfile(1:end-4) '.mat CorrMat ROI file atlas percent0']) %save all Information as .mat
         writematrix(CorrMat, [outputdir filesep  'R2Rconn150_' tempfile(1:end-4) '.csv']) %save Correlation Matrix as .csv 
            
         clear CorrMat i tem* AvgMat AllVolume n* fwdata
    end
end
