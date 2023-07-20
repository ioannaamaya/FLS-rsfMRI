# FLS-rsfMRI

Analysis steps
A. Convert DICOM files to BIDS format including nifti image files and json files. No script provided for this step.  
B. Run the pre-processing pipeline (B0_preprocessing_batch_BIDS)
C.i) Compute functional connectivity (C0_connectivity_batch_BIDS) - ROI-to-ROI FC is step 2 
C.ii) Take average of pre and post scan matrices (average_pre_post)
D. Optional: descriptive analysis of framewise displacement (D1_analyse_FWD)
E. Subtract (contrast between control and experimantal) and average matrices (E0_contrast_and_visualize_batch)
