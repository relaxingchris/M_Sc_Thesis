function [ mass_matrix_trans, stiff_matrix_trans ] = ...
    transform_system(mass_matrix, stiff_matrix, transform_matrix)

% Transformation of mass and stiffness matrix with a given transformation
% matrix.
%
%   Input:
%           Mass matrix             mass_matrix         [kg]
%
%           Stiffness matrix        stiff_matrix        [N/m]
%
%           Transformation matrix   transform_matrix    [-]
%
%   Output:
%           Transformed mass        mass_matrix_trans   [kg]
%           matrix
%
%           Transformed stiffness   stiff_matrix_trans  [N/M]
%           matrix 


mass_matrix_trans = transform_matrix'*mass_matrix*transform_matrix;
stiff_matrix_trans = transform_matrix'*stiff_matrix*transform_matrix;

end

