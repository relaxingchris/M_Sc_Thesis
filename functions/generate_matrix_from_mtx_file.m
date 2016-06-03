function [ output_matrix ] = generate_matrix_from_mtx_file( file_name )

% function that reads a matrix from a .mtx-file and returns
% symmetric matrix containing the values from the file.
%
%   Input:
%           File name               file_name           ['string']
%
%   Output:
%           Symmetric Matrix        output_matrix
%           containing read values


% load file
file_data = dlmread(file_name);

% find rank of matrix from file
rank_matrix = max(max(file_data(:,1:2)));

% generate zero matrix with that rank
output_matrix = zeros(rank_matrix);

% fill matrix
for i=1:size(file_data,1)
    
    % get index where to store the current value
    dof_1 = file_data(i,1);
    dof_2 = file_data(i,2);
    
    output_matrix(dof_1, dof_2) = file_data(i,3);
end

end

