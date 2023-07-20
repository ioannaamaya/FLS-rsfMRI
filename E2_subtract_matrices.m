function E2_subtract_matrices (sub_dir, conn_matrix_pcb, conn_matrix_exp, exp_file_names2, SJs, subject, names)


    for matrix = 1:numel(conn_matrix_pcb)
        A = [conn_matrix_pcb(matrix).folder filesep conn_matrix_pcb(matrix).name];
        matrix_A = load(A);
        if size(conn_matrix_exp) > 0 %if the subject was not excluded in one of the conditions 
        B = [conn_matrix_exp(matrix).folder filesep conn_matrix_exp(matrix).name];
        matrix_B = load(B);
        matrix_C = matrix_B.CorrMat - matrix_A.CorrMat; %subtract placebo from experimental
        ROI = matrix_A.ROI;

        %look for file with NaN values replaced 
        perc_0_file = dir(fullfile(sub_dir, '**', [SJs{subject}, '_average_', exp_file_names2, '_perc0.mat'])); 

        if size(perc_0_file) > 0
            temp = [perc_0_file.folder filesep perc_0_file.name];
            perc_0_matrix = load (temp);
            perc_0_matrix = perc_0_matrix.perc0;
            locations = find(isnan(perc_0_matrix));
            %locations = [x, y];
            matrix_C(locations, :) = NaN; %replace with NaN if zeros > 0.5
        end 

        new_name =  (extractAfter(B, 'R2Rconn170_'));
        new_name =  [(extractBefore(new_name, '.mat'))];
        new_name = [names new_name]; %change name depending on condition
        eval(['save ' sub_dir filesep new_name '.mat matrix_C ROI']);
        writematrix(matrix_C, [sub_dir filesep new_name '.csv']) %save as csv file
        end 
    end    
end


