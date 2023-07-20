function E3_average_matrices (all_matrices, dataset_dir)


all_sub_matrices = zeros((numel(all_matrices)), 150, 150); %initialize big matrix.
for matrix = 1:numel(all_matrices)
    temp = [all_matrices(matrix).folder filesep all_matrices(matrix).name];
    matrix_sub = load(temp);
    all_sub_matrices(matrix, :, :) = matrix_sub.matrix_C;  
    clear matrix_sub;
end
average_matrix = mean(all_sub_matrices, 'omitnan');
average_matrix = squeeze(average_matrix);
std_matrix = std(all_sub_matrices, 'omitnan');
std_matrix = squeeze(std_matrix);
matrix_name = [extractBefore(all_matrices(1).name, "m0.4wrasub")];
matrix_name = ['average_' matrix_name];

%save the mean matrix file as .mat and .csv and one .mat file of all_sub_matrices
eval(['save ' dataset_dir filesep matrix_name '.mat average_matrix']) %save average matrix as .mat file
writematrix(average_matrix, [dataset_dir filesep matrix_name '.csv']) %save as csv 

big_matrix_name = ['all_matrices_test' matrix_name];
    
eval(['save ' dataset_dir filesep big_matrix_name '.mat all_sub_matrices'])
clear sum_matrices;

