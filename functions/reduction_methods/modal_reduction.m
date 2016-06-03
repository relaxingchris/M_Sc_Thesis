function [ ss_sys, M_red, D_damp_red, K_red, T_complete ] = modal_reduction(...
    set_of_modes_kept, M, D_damp, K)

% Function that reduces the input system using Component Mode Synthesis (CMS)
%
%   Input:
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

% compute eigenvectors of full model
[eig_vec_full, ~] = eig_vec_mass_norm_and_eig_val_hz(M, K);

% truncate to kept modes
T_modal = eig_vec_full(:,set_of_modes_kept);

% compute complete transformation matrix
T_complete = T_modal;

% reduce system matrices
M_red = T_complete'*M*T_complete;
D_damp_red = T_complete'*D_damp*T_complete;
K_red = T_complete'*K*T_complete;

% set up explicit state space model
ss_sys = exp_ss_model(M_red, D_damp_red, K_red);

end

