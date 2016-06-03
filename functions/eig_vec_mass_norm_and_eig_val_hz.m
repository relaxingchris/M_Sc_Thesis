function [eig_vec_mass_norm, eig_val_hz ] = ...
    eig_vec_mass_norm_and_eig_val_hz(mass_matrix, stiff_matrix)

% function that returns the mass normalized eigenvectors
% and eigenvalues sorted from lowest to highest
%
%   Input:
%           Stiffness matrix        stiff_matrix        [N/m]
%
%           Mass matrix             mass_matrix         [kg]
%
%   Output:
%           Matrix containing all   eig_vec_mass_norm   [-]
%           mass normalized
%           eigenvectors
%
%           Vector with eigenvalues eig_val_hz          [Hz]
%           sorted from lowest to
%           highest


% eigenvalues and eigenvectors of general problem

[eig_vec, eig_val] = eig(stiff_matrix, mass_matrix, 'qz');

% sort eigenvalues from lowest to highest
% eigenvectors are sorted accordingly

eig_val_array = diag(eig_val);

[eig_val_sort, ind] = sort(eig_val_array);
eig_vec_sort = eig_vec(:,ind);

% eigenvalues in Hz

eig_val_hz = real(sqrt(eig_val_sort)/(2*pi));

% mass normalized eigenvectors

eig_vec_mass_norm = eig_vec_sort/sqrt(eig_vec_sort'*mass_matrix*...
    eig_vec_sort);


end

