function [ ss_sys, M_red, D_damp_red, K_red, T_complete ] = dynamic_reduction(...
    frequ_dynamic_reduction, M, D_damp, K, master_set)

% Function that reduces the input system using dynamic reduction.
%
%   Input:
%           Frequency of the        frequ_dynamic_reduction [rad]
%           dynamic reduction
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

% get dymension of system matrices
dim_of_problem = size(M,2);

% get slave DOF set
slave_set = 1:dim_of_problem;
slave_set(ismember(slave_set,master_set)) = [];

% rearrange mass and stiffness matrix according to master and slave
% DOF set
%M_mm = M(master_set,master_set);
%M_ms = M(master_set,slave_set);
M_sm = M(slave_set,master_set);
M_ss = M(slave_set,slave_set);

%M_sort = [M_mm M_ms;M_sm M_ss];

%K_mm = K(master_set,master_set);
%K_ms = K(master_set,slave_set);
K_sm = K(slave_set,master_set);
K_ss = K(slave_set,slave_set);

%K_sort = [K_mm K_ms;K_sm K_ss];

% compute the dynamic part of transformation matrix
B_ss = -M_ss*power(frequ_dynamic_reduction,2)+K_ss;
B_sm = -M_sm*power(frequ_dynamic_reduction,2)+K_sm;

T_dyn = [eye(size(master_set,2));-inv(B_ss)*B_sm];

% compute complete transformation matrix
T_complete = T_dyn;

% reduce system matrices
M_red = T_complete'*M*T_complete;
D_damp_red = T_complete'*D_damp*T_complete;
K_red = T_complete'*K*T_complete;

% set up explicit state space model
ss_sys = exp_ss_model(M_red, D_damp_red, K_red);
end

