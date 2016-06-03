function [eig_vec_sort, eig_val_hz ] = ...
    eig_vec_and_eig_val_hz(mass_matrix, stiff_matrix)

% function that returns the eigenvectors (NOT MASS NORMALIZED)
% and eigenvalues sorted from lowest to highest
%
%   Input:
%           Stiffness matrix        stiff_matrix        [N/m]
%
%           Mass matrix             mass_matrix         [kg]
%
%   Output:
%           Matrix containing all   eig_vec             [-]
%           eigenvectors
%
%           Vector with eigenvalues eig_val_hz          [Hz]
%           sorted from lowest to
%           highest


% eigenvalues and eigenvectors of general problem
% 'qz' algorithm better suited to solve badly conditioned matrices
[eig_vec, eig_val] = eig(stiff_matrix, mass_matrix, 'qz');

% sort eigenvalues from lowest to highest
% eigenvectors are sorted accordingly

eig_val_array = diag(eig_val);

[eig_val_sort, ind] = sort(eig_val_array);
eig_vec_sort = eig_vec(:,ind);

% eigenvalues in Hz

eig_val_hz = real(sqrt(eig_val_sort)/(2*pi));

end

