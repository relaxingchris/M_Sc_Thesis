function [ exp_ss_sys ] = exp_ss_model( M, D_damp, K )

% Function that returns an explicit MIMO state space system
%
%   Input:
%           Mass matrix             M                   [kg]
%
%           Damping matrix          D_damp              [Ns/m]
%
%           Stiffness matrix        K                   [N/m]
%
%   Output:
%           State space system      exp_ss_sys          [state space system]
%           of reduced problem

% size of problem 
ndof = size(M,2);

A = [zeros(ndof) eye(ndof);-M\K -M\D_damp]; %M\K equals inv(M)*K
B = [zeros(ndof);inv(M)];
C = [eye(ndof) zeros(ndof)];
D = zeros(ndof);

exp_ss_sys = ss(A,B,C,D);

end

