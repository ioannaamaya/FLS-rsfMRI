function E1_perc0_analysis (sub_dir, placebo_matrices, experimental_matrices, name)

   
    for matrix = 1:numel(placebo_matrices)
        temp = [placebo_matrices(matrix).folder filesep placebo_matrices(matrix).name];
        matrix_sub = load (temp);
        perc0_placebo = matrix_sub.Percent0; 
        clear matrix_sub;
        y = find(perc0_placebo>0.5); 
    end 
    if size(y) > 0 
        for nan = 1:size(y)
            perc0_placebo(y(nan)) = NaN; %replace values >0.5 with Nan 
        end 
        ROI = unique(y);
        counts = histc(y(:), ROI);
    end 

    for matrix = 1:numel(experimental_matrices)
        temp = [experimental_matrices(matrix).folder filesep experimental_matrices(matrix).name];
        matrix_sub = load (temp);
        perc0_experimental = matrix_sub.percent0; 
        clear matrix_sub;
        y2 = find(perc0_experimental>0.5); 
    end 
    if size(y2) > 0 
        for nan = 1:size(y2)
            perc0_placebo(y2(nan)) = NaN; %replace values >0.5 with Nan 
        end 
        ROI2 = unique(y2);
        counts2 = histc(y2(:), ROI2);
    end 
    perc0 = perc0_experimental + perc0_placebo; %voxels missing in either experimental or placebo will be coded as NaN. 
    eval(['save ' sub_dir filesep name '.mat perc0']);

