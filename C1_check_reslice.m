function C1_check_reslice(rest_path, atlas)

V = spm_vol(atlas);
U = spm_vol(rest_path);
U_dim = num2str(U(1).dim);
if size(U) == 1
    if V.dim ~= U.dim
        s = struct('mask',{false},'mean',{false},'interp',{0},'which',{2},'wrap',{[1 1 0]},'prefix',{'resliced_'});
        spm_reslice({rest_path; atlas}, s);
    end
else
    if V.dim ~= U(1).dim
        s = struct('mask',{false},'mean',{false},'interp',{0},'which',{2},'wrap',{[1 1 0]},'prefix',{'resliced_'});
        spm_reslice({rest_path; atlas}, s);
    end
end