function [ ss_sys, M_red, D_damp_red, K_red, T_complete ] = CMS_reduction(...
    type_of_cms, set_of_modes_kept, M, D_damp, K, master_set)

% Function that reduces the input system using Component Mode Synthesis (CMS)
%
%   Input:
%           CMS type                type_of_cms         ['string']
%
%           Set containing numbers  set_of_modes_kept   [-]
%           of eigenmodes to be
%           kept
%
%           Mass matrix             M                   [kg]
%
%           Damping matrix          D_damp              [Ns/m]
%
%           Stiffness matrix        K                   [N/m]
%
%           Set containing master   master_set          [-]
%           DOF
%
%   Output:
%           State space system      ss_sys              [state space system]
%           of reduced problem
%
%           Reduced mass matrix     M_red               [kg]
%
%           Reduced damping matrix  D_red               [Ns/m]
%
%           Reduced stiffness       K_red               [N/m]
%           matrix
%
%           Transformation matrix   T                   [-]
%           from reduced back to
%           physical space

if (strcmpi(type_of_cms,'cb')) %CRAIG-BEMPTON REDUCTION
    
    % get dymension of system matrices
    dim_of_problem = size(M,2);
    
    % get slave DOF set
    slave_set = 1:dim_of_problem;
    slave_set(ismember(slave_set,master_set)) = [];
    
    % rearrange mass and stiffness matrix according to master and slave
    % DOF set
    %M_mm = M(master_set,master_set);
    %M_ms = M(master_set,slave_set);
    %M_sm = M(slave_set,master_set);
    M_ss = M(slave_set,slave_set);
    
    %M_sort = [M_mm M_ms;M_sm M_ss];
    
    %K_mm = K(master_set,master_set);
    %K_ms = K(master_set,slave_set);
    K_sm = K(slave_set,master_set);
    K_ss = K(slave_set,slave_set);
    
    %K_sort = [K_mm K_ms;K_sm K_ss];
    
    % compute the static part of transformation matrix
    T_stat = [eye(size(master_set,2));
        -inv(K_ss)*K_sm];
    
    % compute eigenvectors of model corresponding to slave/slave set
    [eig_vec_cb, eig_val_cb] = eig_vec_mass_norm_and_eig_val_hz(M_ss, K_ss);
    
    % truncate to kept modes
    phi_cb = eig_vec_cb(:,set_of_modes_kept);
    fprintf('-----------------------------------------\n')
    fprintf('Eigenfrequency of highest mode in CB set:\n')
    fprintf('%i Hz\n', eig_val_cb(set_of_modes_kept(end)))
    fprintf('-----------------------------------------\n')
    
    % set up CMS part of transformation matrix
    T_cms = [zeros(size(master_set,2),size(set_of_modes_kept,2));
        phi_cb];
    
    % compute complete transformation matrix
    T_complete = [T_stat T_cms];
    
    % reduce system matrices
    M_red = T_complete'*M*T_complete;
    D_damp_red = T_complete'*D_damp*T_complete;
    K_red = T_complete'*K*T_complete;
    
    % set up explicit state space model
    ss_sys = exp_ss_model(M_red, D_damp_red, K_red);
end

end

