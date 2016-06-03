function [ output_vec ] = nat_frequ_diff_criterion( eig_val_full_Hz,...
    eig_val_red_Hz, array_master_DOF_or_kept_modes)

% function that computes the modified modal assurance criterion of
% two eigenvector matrices with mass normalization.
% After computation, the MAC will be plotted.
%
%   Input:
%           Matrix with full        eig_vec_full        [-]
%           eigenvectors
%
%           Matrix with reduced     eig_vec_red         [-]
%           eigenvectors
%
%           Mass matrix (not        mass_matrix         [kg]
%           reduced) 
%
%           Array with DOF number   array_Master_DOF    [-]
%           that are definded       _or_kept_modes
%           as Master nodes
%           or modes that are
%           kept
%
%           'real'/'complex'        type                string  
%           eigenvectors
%
%   Output:
%           Matrix containing the   output_matrix       [-]
%           modal assurance values

% pre allocate return vector
output_vec = zeros(size(eig_val_red_Hz,1),1);

% compute each entry of output vector
for i=1:size(eig_val_red_Hz,1)
    
    master_DOF = array_master_DOF_or_kept_modes(i);
    
    if abs(eig_val_full_Hz(master_DOF)) < 1e-3
        output_vec(i) = NaN;
    else
        output_vec(i) = abs(eig_val_red_Hz(i) - eig_val_full_Hz(master_DOF)) / ...
            eig_val_full_Hz(master_DOF);
    end
    
end

% plot result

figure
h = plot(1:size(eig_val_red_Hz,1),output_vec*100);
title('Natural Frequency Difference')
ylabel('rel. error [%]')

end


