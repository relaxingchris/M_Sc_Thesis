function [ imp_ss_sys ] = imp_ss_model( M, D_damp, K )

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
%           State space system      imp_ss_sys          [state space system]
%           of reduced problem

% size of problem 
ndof = size(M,2);

A = [zeros(ndof) eye(ndof);-K -D_damp]; %M\K equals inv(M)*K
B = [zeros(ndof);eye(ndof)];
C = [eye(ndof) zeros(ndof)];
D = zeros(ndof);
E = [eye(ndof) zeros(ndof);zeros(ndof) M];

imp_ss_sys = dss(A,B,C,D,E);

end

