function E4_plot_conn_matrix (dataset_dir, name, suffix)

plot_range = [1:11 13:22 24:36 38:47 49:50 121:150]; % FOR HYBRID ATLAS


ROIs = textread('ROI_lat.txt', '%s', 'delimiter','\n'); % for hybrid atlas

colour_scale =  [0,0.05,1;
        0,0.01,1;
        0,0.02,1;
        0,0.25,1;
        0,0.3,1;
        0,0.4,1;
        0,0.50,1;
        0,0.6,1;
        0,0.7,1;
        0,0.75,1;
        0,0.8,1;
        0,0.9,1;
        0,1,1;
        1,1,1;
        1,1,0;
        1,0.9,0;
        1,0.8,0;
        1,0.75,0;
        1,0.7,0;
        1,0.6,0;
        1,0.5,0;
        1,0.4,0;
        1,0.3,0;
        1,0.25,0;
        1,0.2,0;
        1,0.1,0;
        1,0.050,0];


average_matrix = dir(fullfile(dataset_dir, '**', ['average_', name, suffix]));
temp = [average_matrix.folder filesep average_matrix.name];
average_matrix2 = load (temp);
average_matrix3 = average_matrix2.average_matrix(plot_range,plot_range);


xlswrite('mtx_10Hz.xlsx', average_matrix3)

name = [extractBefore(average_matrix.name, ".mat")];
name = ['cut_ROI' name]; 
%plot conn matrix 
clims = [-0.4 0.4];
image = imagesc(average_matrix3, clims);
colormap(colour_scale);
cb = colorbar;
%save image 
cd ([dataset_dir filesep 'figures'])
saveas(image, name, 'jpg');





 thal_idx = 121:150;
 cort_idx = [1:23];

ROI_labels = string(ROIs); 

clims = [-0.4 0.4];

 set(gca, 'XTick', [1:length(plot_range)], 'XTickLabel', ROI_labels, 'FontSize', 7) % or put cortical_aal3 for hybrid
 set(gca, 'YTick', [1:length(plot_range)], 'YTickLabel', ROI_labels, 'FontSize', 7) % or put thalamus_aal3 for hybrid
