function compute_conn_map(nii_fmri_in_path, nii_seed_mask_in_path, nii_searchlight_mask_in_path, output_prefix)
%
%  compute_conn_map(nii_fmri_in_path, nii_seed_mask_in_path, nii_searchlight_mask_in_path, output_prefix)
%
%  Computes connectivity (i.e. correlation) maps for each voxel in a given searchlight with respect to a given seed.
%
%  Requires MATLAB R2017b or later.
%
%  Required MATLAB library:
%
%  Tools for NIfTI and ANALYZE image (https://www.mathworks.com/matlabcentral/fileexchange/8797-tools-for-nifti-and-analyze-image).
%
%  Input MATLAB variables:
%
%    - nii_fmri_in_path               : (Required) Path to the Nifti file of the fMRI timeseries, as a string or character array.
%    - nii_seed_mask_in_path          : (Required) Path to the Nifti file of the seed mask, as a string or character array. Voxels with a value
%                                       greater than zero are considered as being inside the mask.
%    - nii_searchlight_mask_in_path   : (Required) Path to the Nifti file of the searchlight mask, as a string or character array. Voxels with
%                                       a value greater than zero are considered as being inside the mask.
%    - output_prefix                  : (Required) Prefix to prepend on the name of the output Nifti files, as a string or character array. A
%                                       full or relative path can be included with the prefix (e.g. /path/to/file/prefix).
%
%
%  Resulting Nifti files:
%
%    - {output_prefix}conn_map.nii.gz : Nifti file containing the Z-transformed connectivity map.
%
%                       
%  (c) 2022-2023. Written by Dr. Stefan Pszczolkowski (mszspp@nottingham.ac.uk) 
%
    if nargin < 4
        error('Not enough input arguments.');
    end
    
    output_prefix = strrep(output_prefix, '\', '/');
    prefix_path = fileparts(output_prefix);
    if ~isempty(prefix_path) && ~exist(prefix_path, 'dir')
        error('Output folder ''%s'' does not exist.', prefix_path); 
    end
    
    if ~exist(nii_fmri_in_path, 'file')
       error('File ''%s'' does not exist.', nii_fmri_in_path); 
    end
    
    if ~exist(nii_seed_mask_in_path, 'file')
       error('File ''%s'' does not exist.', nii_seed_mask_in_path); 
    end
    
    if ~exist(nii_searchlight_mask_in_path, 'file')
       error('File ''%s'' does not exist.', nii_searchlight_mask_in_path); 
    end
    
    % Make sure string-type input are (converted to) character arrays
    output_prefix = convertStringsToChars(output_prefix);
    nii_fmri_in_path = convertStringsToChars(nii_fmri_in_path);
    nii_seed_mask_in_path = convertStringsToChars(nii_seed_mask_in_path);
    nii_searchlight_mask_in_path = convertStringsToChars(nii_searchlight_mask_in_path);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    fprintf('Loading images ... ');
    nii_fmri = load_untouch_nii(nii_fmri_in_path);
    nii_mask_seed = load_untouch_nii(nii_seed_mask_in_path);
    nii_mask_searchlight = load_untouch_nii(nii_searchlight_mask_in_path);
    fprintf('done\n');
    
    [timeseries_size_x, timeseries_size_y, timeseries_size_z, timeseries_size_t] = size(nii_fmri.img);
    
    if ~isequal(size(nii_mask_seed.img), [timeseries_size_x, timeseries_size_y, timeseries_size_z])
       error('The size of the seed mask image does not match the volume size of the fMRI timeseries.'); 
    end
    
    if ~isequal(size(nii_mask_searchlight.img), [timeseries_size_x, timeseries_size_y, timeseries_size_z])
       error('The size of the searchlight mask image does not match the volume size of the fMRI timeseries.'); 
    end

    fprintf('Initialising timeseries ... ');
    volume_size = [timeseries_size_x, timeseries_size_y, timeseries_size_z];
    volume_length = prod(volume_size);

    timepoint_data = reshape(nii_fmri.img, [volume_length, timeseries_size_t]);
    score_map = zeros(volume_length, 1);
    
    mask_brain = (var(timepoint_data, [], 2) > 0);
    
    mask_seed = (nii_mask_seed.img(:) > 0) & mask_brain;
    mask_searchlight = (nii_mask_searchlight.img(:) > 0) & mask_brain & (~mask_seed);
    
    X = mean(timepoint_data(mask_seed,:))';
    fprintf('done\n');
    
    fprintf('Computing correlations ... ');
    searchlight_indices = find(mask_searchlight)';
    
    for index = searchlight_indices   
        Y = timepoint_data(index,:)';

        score_map(index) = corr(X, Y);
    end
    fprintf('done\n');

    fprintf('Saving image ... ');
    
    nii_out = nii_mask_seed;
    nii_out.img = reshape(score_map, [timeseries_size_x, timeseries_size_y, timeseries_size_z]);
    nii_out.hdr.dime.bitpix = 32;
    nii_out.hdr.dime.datatype = 16;
    nii_out.hdr.dime.glmin = min(nii_out.img(:));
    nii_out.hdr.dime.glmax = max(nii_out.img(:));
    save_untouch_nii(nii_out, [output_prefix 'conn_map.nii.gz']);
    
    fprintf('done\n');
end

