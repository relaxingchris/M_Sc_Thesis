function [ output_matrix ] = MAC_mass_mod( eig_vec_1, eig_vec_2, ...
    mass_matrix, array_master_DOF_or_kept_modes, type)

% function that computes the modified modal assurance criterion of
% two eigenvector matrices with mass normalization.
% After computation, the MAC will be plotted.
%
%   Input:
%           Matrix1 with            eig_vec_1           [-]
%           eigenvectors
%
%           Matrix2 with            eig_vec_2           [-]
%           eigenvectors
%
%           Mass matrix (not        mass_matrix         [kg]
%           reduced) 
%
%           Array with DOF number   array_master_DOF_   [-]
%           that are definded       or_kept_modes
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

% pre allocate return matrix
output_matrix = zeros(size(eig_vec_2,2));

% if only real parts should be used
if strcmpi(type,'real')
    fprintf(['Only real parts of eigenvectors are used for MAC ', ...
    'computation!\n']);
    eig_vec_1 = real(eig_vec_1);
    eig_vec_2 = real(eig_vec_2);
elseif strcmpi(type,'complex')
    fprintf(['Real and complex parts of eigenvectors are used for MAC ', ...
    'computation!\n']);
else
    msg = 'Input string not equal to ''real'' or ''complex''.';
    error(msg)
end

% calculate all single MAC values with loop
for i=1:size(eig_vec_2,2)
    master_DOF = array_master_DOF_or_kept_modes(i);
    for j=1:size(eig_vec_2,2)
        output_matrix(i,j)=single_mac(eig_vec_1(:,master_DOF),...
            eig_vec_2(:,j), mass_matrix);
    end
end

%% plot mac matrix with color of each bar according to "height"
figure

h = bar3(output_matrix);
view(-70,50);

for i = 1:size(output_matrix,2)
    cdata = get(h(i),'cdata');
    k = 1;
    for j = 0:6:(6*size(output_matrix,1)-6)
        cdata(j+1:j+6,:) = output_matrix(k,i);
        k = k+1;
    end
    set(h(i),'cdata',cdata);
end

% change colormap
colormap(jet);

% enable colorbar
colorbar

title('mass normalized MAC')

end

function retrn = single_mac(vec_1,vec_2,M)
% This function calculates the modal modified MAC between vector
% phi1 and phi2
retrn = (vec_1'*M*vec_2)^2/((vec_1'*M*vec_1)*(vec_2'*M*vec_2));
end
