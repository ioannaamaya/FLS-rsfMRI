clc
clear all

%data source directories
src_dirs = {'F:\BIDS_Flicker\derivatives'};
for dataset = 1:size(src_dirs, 1)
    src_dir = src_dirs{dataset};

all_files = dir(fullfile(src_dir, 'sub*')); % looking for subject's folders
SJs = {all_files.name}; % creating cell array with subjects' IDs 

for sj = 1:numel(SJs)
    sub_dir = [all_files(sj).folder filesep all_files(sj).name];
    fwd_files = dir(fullfile(sub_dir, '**', ['*bold_FWDstat.mat'])); 
    rp_files = dir(fullfile(sub_dir, '**', ['rp_asub', '*.txt']));
    fwd_names = {fwd_files.name};
    tasks = unique(extractAfter(fwd_names, "task-"));
    tasks = (extractBefore(tasks, "_"))'; %here we extract the task names directly from the file names.
        for file = 1:numel(fwd_files)
            current_fwd = [fwd_files(file).folder filesep fwd_files(file).name];
            current_rp = [rp_files(file).folder filesep rp_files(file).name];
            motion_stats = load(current_rp); 
            load((current_fwd),'fwd');
            mFWD(sj,file)=mean(fwd((segment_start(sj)):(segment_end(sj)))); %save mean fwd for each subject and each session (#file = #ses)
            sFWD(sj,file)=std(fwd((segment_start(sj)):(segment_end(sj))));
            maxFWD(sj,file)=max(fwd);
            accFWD(sj, file)=sum(fwd);
            perc1(sj,file)=(sum(fwd>0.1)/length(fwd))*100;
            perc2(sj,file)=(sum(fwd>0.2)/length(fwd))*100;
            perc3(sj,file)=(sum(fwd>0.3)/length(fwd))*100;
            perc4(sj,file)=(sum(fwd>0.4)/length(fwd))*100;
            perc5(sj,file)=(sum(fwd>0.5)/length(fwd))*100;
            acc_x(sj,file)=sum(diff(motion_stats(:, 1)));
            acc_y(sj,file)=sum(diff(motion_stats(:, 2)));
            acc_z(sj,file)=sum(diff(motion_stats(:, 3)));
            acc_rot_x(sj,file)=sum(diff(motion_stats(:, 4)));
            acc_rot_y(sj,file)=sum(diff(motion_stats(:, 5)));
            acc_rot_z(sj,file)=sum(diff(motion_stats(:, 6)));
            acc_xyz(sj,file)=sum(sum(diff(motion_stats(:, 1:3))));
            acc_rotation(sj,file)=sum(sum(diff(motion_stats(:, 4:6))));
           % motion_3mm_x(sj, file) = (sum(diff(motion_stats(:,1)) > 0.3)) / (length(motion_stats))* 100;
            %motion_3mm_y(sj, file) = (sum(diff(motion_stats(:,2)) > 0.3)) / (length(motion_stats)) *100;
            %motion_3mm_z(sj, file) = (sum(diff(motion_stats(:, 3)) > 0.3)) / (length(motion_stats))*100;
            %motion_3mm_rx(sj, file) = (sum(diff(motion_stats(:, 4)) > 0.3)) / (length(motion_stats))*100;
            %motion_3mm_ry(sj, file) = (sum(diff(motion_stats(:, 5)) > 0.3)) / (length(motion_stats)) * 100;
            %motion_3mm_rz(sj, file) = (sum(diff(motion_stats(:, 6)) > 0.3)) / (length(motion_stats)) * 100;

            clear fwd
        end
end

table = [mFWD sFWD];
%table = [mFWD sFWD maxFWD accFWD perc3 perc4 perc5 acc_x acc_y acc_z acc_rot_x acc_rot_y acc_rot_z acc_xyz acc_rotation];
all_FWD = array2table(table);
%all_FWD = array2table(table, 'VariableNames', {'mFWD run1', 'mFWD run2', 'sFWD run1', 'sFWD run2', 'maxFWD run1', 'maxFWD run2', 'accFWD run1', 'accFWD run2', 'perc1 run1', 'perc1 run 2', 'perc2 run 1', 'perc2 run 2', 'perc3 run 1', 'perc3 run2'});
eval(['save ' src_dir filesep 'all_FWD_N2' '.mat all_FWD']);
fig = figure;
plot(mFWD);
legend(tasks);
xlabel('subject number')
ylabel('mean fwd')
cd (src_dir);
saveas(fig, 'plot_mean_FWD', 'fig');
clear fig
fig = figure;
plot(sFWD);
legend(tasks);
xlabel('subject number')
ylabel('sd fwd')
cd (src_dir);
saveas(fig, 'plot_mean_sWD', 'fig');
writetable(all_FWD, 'allFWD.xlsx')
end

%%
mFWD([1 4],:)=[];
[p,anovatab,stats] = anova1(mFWD);

[hV1 pV1]=ttest(mFWD(:,1),mFWD(:,2)); % pre vs ganz
[hV2 pV2]=ttest(mFWD(:,2),mFWD(:,3)); % ganz vs post
[hV3 pV3]=ttest(mFWD(:,1),mFWD(:,3)); % pre vs post


[p,anovatab,stats] = anova1(sFWD);

[hV1 pV1]=ttest(sFWD(:,1),sFWD(:,2)); % pre vs ganz
[hV2 pV2]=ttest(sFWD(:,2),sFWD(:,3)); % ganz vs post
[hV3 pV3]=ttest(sFWD(:,1),sFWD(:,3)); % pre vs post
