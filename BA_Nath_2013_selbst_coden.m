%Selbst Coden von Bachelorarbeit von Nath, 2013
%
%

clear
clc

M = diag([1 2 2 2 1]);
K = [100 -100 0 0 0;
    -100 200 -100 0 0;
    0 -100 200 -100 0;
    0 0 -100 200 -100;
    0 0 0 -100 100];

% eigenvalues and eigenvectors of general problem

[eigVec, eigVal] = eig(K,M);

% eigenvalues in Hz

eigVecArray = diag(eigVal);
eigValHz = sqrt(eigVecArray)/(2*pi);

% mass normalized eigenvectors

eigVecMassNorm = eigVec/sqrt(eigVec'*M*eigVec);

% calculation of modal mass and stiffness matrix

M_mod = eigVecMassNorm'*M*eigVecMassNorm;
K_mod = eigVecMassNorm'*K*eigVecMassNorm;

%% modal reduction

% only contain first 3 modes
T_modal_redu = eigVecMassNorm(:,1:3);

M_modal_redu = T_modal_redu'*M*T_modal_redu;
K_modal_redu = T_modal_redu'*K*T_modal_redu;

%% Guyan with 2 external nodes (DOF 1 and 5)

% reorganize mass and stiffness matrix
M_guy2_hh = M([1,5],[1,5]);
M_guy2_hn = M([1,5],[2:4]);
M_guy2_nh = M([2:4],[1,5]);
M_guy2_nn = M([2:4],[2:4]);

M_guy2_sort = [M_guy2_hh M_guy2_hn;
    M_guy2_nh M_guy2_nn];

K_guy2_hh = K([1,5],[1,5]);
K_guy2_hn = K([1,5],[2:4]);
K_guy2_nh = K([2:4],[1,5]);
K_guy2_nn = K([2:4],[2:4]);

K_guy2_sort = [K_guy2_hh K_guy2_hn;
    K_guy2_nh K_guy2_nn];

% static transformation matrix

T_s_guy2 = -inv(K_guy2_nn)*K_guy2_nh;

% full transformation matrix

T_guy2 = [eye(2);T_s_guy2];

%% Guyan with 3 external nodes (DOF 1, 3 and 5)

% reorganize mass and stiffness matrix
M_guy3_hh = M([1,3,5],[1,3,5]);
M_guy3_hn = M([1,3,5],[2,4]);
M_guy3_nh = M([2,4],[1,3,5]);
M_guy3_nn = M([2,4],[2,4]);

M_guy3_sort = [M_guy3_hh M_guy3_hn;
    M_guy3_nh M_guy3_nn];

K_guy3_hh = K([1,3,5],[1,3,5]);
K_guy3_hn = K([1,3,5],[2,4]);
K_guy3_nh = K([2,4],[1,3,5]);
K_guy3_nn = K([2,4],[2,4]);

K_guy3_sort = [K_guy3_hh K_guy3_hn;
    K_guy3_nh K_guy3_nn];

% static transformation matrix

T_s_guy3 = -inv(K_guy3_nn)*K_guy3_nh;

% full transformation matrix

T_guy3 = [eye(3);T_s_guy3];

%%  Craig Bempton Method with 2 external nodes (DOF 1 and 5) and
%   1 additional mode

% reorganize mass and stiffness matrix
M_cb1_hh = M([1,5],[1,5]);
M_cb1_hn = M([1,5],[2:4]);
M_cb1_nh = M([2:4],[1,5]);
M_cb1_nn = M([2:4],[2:4]);

M_cb1_sort = [M_cb1_hh M_cb1_hn;
    M_cb1_nh M_cb1_nn];

K_cb1_hh = K([1,5],[1,5]);
K_cb1_hn = K([1,5],[2:4]);
K_cb1_nh = K([2:4],[1,5]);
K_cb1_nn = K([2:4],[2:4]);

K_cb1_sort = [K_cb1_hh K_cb1_hn;
    K_cb1_nh K_cb1_nn];

% static transformation matrix

T_s_cb1 = -inv(K_cb1_nn)*K_cb1_nh;

% dynamic transformation matrix

[eig_vec_cb1_CB, eig_val_cb1_CB] = eig(K_cb1_nn,M_cb1_nn);

% eigenvalues in Hz

% eig_val_array_cb1 = diag(eig_val_cb1);
% eig_val_cb1_Hz = sqrt(eig_val_array_cb1)/(2*pi);

% mass normalized eigenvectors

eig_vec_mass_norm_cb1 = eig_vec_cb1_CB/sqrt(eig_vec_cb1_CB'*M_cb1_nn*eig_vec_cb1_CB);

% full transformation matrix

T_cb1 = [eye(2) zeros(2,1);
    T_s_cb1 eig_vec_mass_norm_cb1(:,1)];


%%  Craig Bempton Method with 2 external nodes (DOF 1 and 5) and
%   2 additional modes

% reorganize mass and stiffness matrix
M_cb2_hh = M([1,5],[1,5]);
M_cb2_hn = M([1,5],[2:4]);
M_cb2_nh = M([2:4],[1,5]);
M_cb2_nn = M([2:4],[2:4]);

M_cb2_sort = [M_cb2_hh M_cb2_hn;
    M_cb2_nh M_cb2_nn];

K_cb2_hh = K([1,5],[1,5]);
K_cb2_hn = K([1,5],[2:4]);
K_cb2_nh = K([2:4],[1,5]);
K_cb2_nn = K([2:4],[2:4]);

K_cb2_sort = [K_cb2_hh K_cb2_hn;
    K_cb2_nh K_cb2_nn];

% static transformation matrix

T_s_cb2 = -inv(K_cb2_nn)*K_cb2_nh;

% dynamic transformation matrix

[eig_vec_cb2_CB, eig_val_cb2_CB] = eig(K_cb2_nn,M_cb2_nn);

% eigenvalues in Hz

% eig_val_array_cb1 = diag(eig_val_cb1);
% eig_val_cb1_Hz = sqrt(eig_val_array_cb1)/(2*pi);

% mass normalized eigenvectors

eig_vec_mass_norm_cb2 = eig_vec_cb2_CB/sqrt(eig_vec_cb2_CB'*...
    M_cb2_nn*eig_vec_cb2_CB);

% full transformation matrix

T_cb2 = [eye(2) zeros(2,2);
    T_s_cb2 eig_vec_mass_norm_cb2(:,1:2)];

%% get eigenfrequencies of reduced models

% ORIGINAL SYSTEM

fprintf('Original System:\n')
eigValHz

% MODAL REDUCTION

% setting up reduced mass and stiffness matrices

[M_mod_red, K_mod_red] = transform_system(M, K, T_modal_redu);

[eig_vec_mod_red, eig_val_mod_red] = eig_vec_mass_norm_and_eig_val_hz(...
    M_mod_red, K_mod_red);

fprintf('Modal reduction System:\n')
eig_val_mod_red


% GUYAN WITH 2 DOF (DOF 1 and 5)

% setting up reduced mass and stiffness matrices

[M_guy2_red, K_guy2_red] = transform_system(M_guy2_sort, K_guy2_sort,...
    T_guy2);

[eig_vec_guy2_red, eig_val_guy2_red] = eig_vec_mass_norm_and_eig_val_hz(...
    M_guy2_red, K_guy2_red);

fprintf('Guyan 2 System:\n')
eig_val_guy2_red


% GUYAN WITH 3 DOF (DOF 1, 3 and 5)

% setting up reduced mass and stiffness matrices

[M_guy3_red, K_guy3_red] = transform_system(M_guy3_sort, K_guy3_sort,...
    T_guy3);

[eig_vec_guy3_red, eig_val_guy3_red] = eig_vec_mass_norm_and_eig_val_hz(...
    M_guy3_red, K_guy3_red);

fprintf('Guyan 3 System:\n')
eig_val_guy3_red


% CB WITH 2 DOF AND 1 MODE (DOF 1, 3 and 5)

% setting up reduced mass and stiffness matrices

[M_cb1_red, K_cb1_red] = transform_system(M_cb1_sort, K_cb1_sort,...
    T_cb1);

[eig_vec_cb1_red, eig_val_cb1_red] = eig_vec_mass_norm_and_eig_val_hz(...
    M_cb1_red, K_cb1_red);

fprintf('CMS with 2 external Nodes and 1 additional mode:\n')
eig_val_cb1_red


% CB WITH 2 DOF AND 2 MODES (DOF 1, 3 and 5)

% setting up reduced mass and stiffness matrices

[M_cb2_red, K_cb2_red] = transform_system(M_cb2_sort, K_cb2_sort,...
    T_cb2);

[eig_vec_cb2_red, eig_val_cb2_red] = eig_vec_mass_norm_and_eig_val_hz(...
    M_cb2_red, K_cb2_red);

fprintf('CMS with 2 external Nodes and 2 additional modes:\n')
eig_val_cb2_red

%% mass modified MAC

% transform into full domain
eig_vec_mod_red_retrans = T_modal_redu*eig_vec_mod_red;
eig_vec_cb2_red_retrans = T_cb2*eig_vec_cb2_red;

MAC_mass_mod(eigVecMassNorm, eig_vec_mod_red_retrans, M,...
    1:size(eig_vec_mod_red_retrans,2), 'real')
MAC_mass_mod(eigVecMassNorm, eig_vec_cb2_red_retrans, M,...
    1:size(eig_vec_cb2_red_retrans,2), 'real')





