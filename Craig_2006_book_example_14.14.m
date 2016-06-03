% clear all
clear
clc

% define symbolic variables
syms E I L rho A
E_val = 2.1e11;
I_val = 0.02*0.003^3/12;
L_val = 0.45/2;
rho_val = 7850;
A_val = 0.02*0.003;

% set up full stiffness and mass matrix
K = E*I/(L^3)*[12+12 -12 6*L-6*L 6*L;
    -12 12 -6*L -6*L;
    6*L-6*L -6*L 4*L^2+4*L^2 2*L^2;
    6*L -6*L 2*L^2 4*L^2];

M = rho*A*L/420*[156+156 54 22*L-22*L -13*L;
    54 156 13*L -22*L;
    22*L-22*L 13*L 4*L^2+4*L^2 -3*L^2;
    -13*L -22*L -3*L^2 4*L^2];

% substitute real values into symbolic variables
% and convert matrix from symbolic to double
K_val = double(subs(K,[E,I,L,rho,A],[E_val,I_val,L_val,rho_val,A_val]));
M_val = double(subs(M,[E,I,L,rho,A],[E_val,I_val,L_val,rho_val,A_val]));


% define master and slave matrices
master_DOF = 3:4;
slave_DOF = 1:rank(M);
slave_DOF(ismember(slave_DOF,master_DOF)) = [];

M_hh = M_val(master_DOF,master_DOF); 
M_hn = M_val(master_DOF,slave_DOF);
M_nh = M_val(slave_DOF,master_DOF);
M_nn = M_val(slave_DOF,slave_DOF);

K_hh = K_val(master_DOF,master_DOF); 
K_hn = K_val(master_DOF,slave_DOF);
K_nh = K_val(slave_DOF,master_DOF);
K_nn = K_val(slave_DOF,slave_DOF);

% compute Guyan transformation matrix
T_guy = [eye(size(master_DOF,2));-inv(K_nn)*K_nh];

% compute reduced stiffness and mass matrix
K_red = T_guy'*K_val*T_guy;
M_red = T_guy'*M_val*T_guy;

% compute eigenvectors and eigenvalues of full problem
[eig_vec_full, eig_val_full] = eig_vec_mass_norm_and_eig_val_hz(...
    M_val, K_val);

% compute eigenvectors and eigenvalues of reduced problem
[eig_vec_red, eig_val_red] = eig_vec_mass_norm_and_eig_val_hz(...
    M_red, K_red);

%% compare to lumped mass system

K_lumped_red = 48*E*I/(7*L^3)*[16 -5;-5 2];
M_lumped_red = rho*A*L/2*[2 0;0 1];

K_lumped_red_val = double(subs(K_lumped_red,[E,I,L,rho,A],...
    [E_val,I_val,L_val,rho_val,A_val]));
M_lumped_red_val = double(subs(M_lumped_red,[E,I,L,rho,A],...
    [E_val,I_val,L_val,rho_val,A_val]));

[eig_vec_lumped_red, eig_val_lumped_red] = ...
    eig_vec_mass_norm_and_eig_val_hz(M_lumped_red_val, K_lumped_red_val);

%% compare to analytical eigenfrequencies

eig_frequ_factor = [1.8751; 4.69409; 7.8539; 10.99557; 14.1372; 17.279];
eig_frequ_analytic = eig_frequ_factor.^2.*sqrt(E_val*I_val/...
    (rho_val*A_val*(2*L_val)^4))/(2*pi);

%% display frequencies of all versions

fprintf('Eigenfrequencies of analytical solution:\n')
%fprintf('%s\n', eig_frequ_analytic(1:4));
eig_frequ_analytic(1:4)'

fprintf('Eigenfrequencies of full discrete solution:\n')
eig_val_full'

fprintf('Eigenfrequencies of reduced discrete solution:\n')
eig_val_red'

fprintf('Eigenfrequencies of reduced lumped mass solution:\n')
eig_val_lumped_red'

%% mass normalized MAC

% transform eigenvector into full space
eig_vec_red_retrans = T_guy*eig_vec_red;

% mass normalized MAC
MAC_full_red = MAC_mass_mod(eig_vec_full, eig_vec_red_retrans, M_val, ...
    master_DOF, 'real')

%% eigenfrequency comparison

NFD_full_red = nat_frequ_diff_criterion( eig_val_full,...
    eig_val_red, master_DOF);

