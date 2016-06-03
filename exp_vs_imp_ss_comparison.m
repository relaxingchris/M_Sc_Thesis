clear all
close all
clc

addpath('/Users/LaKrizz/Documents/MATLAB/Master_Thesis/model_data')

M = generate_matrix_from_mtx_file('Bar_5_elements_fixed_CB_start_MASS2.mtx');
K = generate_matrix_from_mtx_file('Bar_5_elements_fixed_CB_start_STIF2.mtx');

% unfixed problem more detailed
% M = generate_matrix_from_mtx_file('Beam_3_element_rectangle_COORDINATE_MASS2.mtx');
% K = generate_matrix_from_mtx_file('Beam_3_element_rectangle_COORDINATE_STIF2.mtx');

% remove rigid body modes
BC_set = [1:9];
K = impose_BC_on_stiffness_matrix(K, BC_set);

% remove values close to zero
%K(K<1e-5)=0;

% get number of DOF
ndof = rank(M);

% construct damping
C_d = 0.2*M;

%C_d = zeros(ndof);

%imlicit state space representation

E=[eye(ndof) zeros(ndof);zeros(ndof) M];
A_imp=[zeros(ndof) eye(ndof);-K -C_d];
B_imp=[zeros(ndof);eye(ndof)];
C=[eye(ndof) zeros(ndof)];
D=zeros(ndof);
sys_imp=dss(A_imp,B_imp,C,D,E);
sysss_imp = sss(sys_imp)

% explicit state space representation

A = [zeros(ndof) eye(ndof);-M\K -M\C_d]; %M\K equals inv(M)*K
B = [zeros(ndof);inv(M)];
sys_exp = ss(A,B,C,D);
sysss_exp = sss(sys_exp)

% truncate zu SISO
sysss_imp_trunc = truncate(sysss_imp,10,15)
sysss_exp_trunc = truncate(sysss_exp,10,15)