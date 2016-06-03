% CMS
close all
clear
clc

%addpath('/Users/LaKrizz/Documents/MATLAB/Master_Thesis/model_data')

M = generate_matrix_from_mtx_file(...
    '/Users/LaKrizz/Documents/MATLAB/Master_Thesis/model_data/Bar_5_elements_fixed_CB_start_MASS2.mtx');
K = generate_matrix_from_mtx_file(...
    '/Users/LaKrizz/Documents/MATLAB/Master_Thesis/model_data/Bar_5_elements_fixed_CB_start_STIF2.mtx');

% construct damping
%D_damp = 0.2*M;
D_damp=zeros(rank(M));

% main parameters for methods
master_DOF = [1:12]; %for CMS and Dynamic reduction
modes_to_keep_CMS_1 = 1:20;
modes_to_keep_CMS_2 = 1:40;
frequ_dynamic_reduction = 3595; %in rad (ca. 3595.4 is a eigenvalue)

% parameters for FRF comparison
output_DOF = 1;
input_DOF = 1;
sampling_frequ = 3.594e3:.001:3.596e3; %in rad

%% Master-Slave set generation
% set up slave DOF
slave_DOF = 1:rank(M);
slave_DOF(ismember(slave_DOF,master_DOF)) = [];

% sort matrices
M_mm = M(master_DOF,master_DOF);
M_ms = M(master_DOF,slave_DOF);
M_sm = M(slave_DOF,master_DOF);
M_ss = M(slave_DOF,slave_DOF);

M_sort = [M_mm M_ms;
    M_sm M_ss];

K_mm = K(master_DOF,master_DOF);
K_ms = K(master_DOF,slave_DOF);
K_sm = K(slave_DOF,master_DOF);
K_ss = K(slave_DOF,slave_DOF);

K_sort = [K_mm K_ms;
    K_sm K_ss];


%% compute the static part of the transformation matrix
T_stat = [eye(size(master_DOF,2));
        -inv(K_ss)*K_sm];

%% compute the dynamic transformation matrix
B_ss = -M_ss*power(frequ_dynamic_reduction,2)+K_ss;
B_sm = -M_sm*power(frequ_dynamic_reduction,2)+K_sm;
T_dyn = [eye(size(master_DOF,2));
        -inv(B_ss)*B_sm];
    
%% compute the CMS (i.e. CB) transformation matrix
% set up CB modes
[eig_vec_CB, eig_val_CB] = eig_vec_mass_norm_and_eig_val_hz(M_ss, K_ss);

% truncate to kept modes
phi_CB_CMS_1 = eig_vec_CB(:,modes_to_keep_CMS_1);
phi_CB_CMS_2 = eig_vec_CB(:,modes_to_keep_CMS_2);

T_CB_CMS_1 = [zeros(size(master_DOF,2),size(modes_to_keep_CMS_1,2));
    phi_CB_CMS_1];
T_CB_CMS_2 = [zeros(size(master_DOF,2),size(modes_to_keep_CMS_2,2));
    phi_CB_CMS_2]; 

% compute complete transformation matrix
T_CMS_1 = [T_stat T_CB_CMS_1];
T_CMS_2 = [T_stat T_CB_CMS_2];

%% FRF comparison

% full system
ndof_full = rank(M);

A_full = [zeros(ndof_full) eye(ndof_full);-M\K -M\D_damp]; %M\K equals inv(M)*K
B_full = [zeros(ndof_full);inv(M)];
C_full = [eye(ndof_full) zeros(ndof_full)];
D_full = zeros(ndof_full);

sys_full = ss(A_full,B_full,C_full,D_full); %create MIMO state space system
sysss_full = sss(sys_full); %convert into sparse state space system
sysss_full_trunc = truncate(sysss_full,output_DOF,input_DOF); %reduce to SISO

% reduced system dynamic reduction
M_red_dyn = T_dyn'*M*T_dyn;
K_red_dyn = T_dyn'*K*T_dyn;
D_damp_red_dyn = T_dyn'*D_damp*T_dyn;

ndof_red_dyn = rank(M_red_dyn);

A_red_dyn = [zeros(ndof_red_dyn) eye(ndof_red_dyn);-M_red_dyn\...
    K_red_dyn -M_red_dyn\D_damp_red_dyn]; %M\K equals inv(M)*K
B_red_dyn = [zeros(ndof_red_dyn);inv(M_red_dyn)];
C_red_dyn = [eye(ndof_red_dyn) zeros(ndof_red_dyn)];
D_red_dyn = zeros(ndof_red_dyn);

sys_red_dyn = ss(A_red_dyn,B_red_dyn,C_red_dyn,...
    D_red_dyn); %create MIMO state space system
sysss_red_dyn = sss(sys_red_dyn); %convert into sparse state space system
sysss_red_dyn_trunc = truncate(sysss_red_dyn,...
    output_DOF,input_DOF); %reduce to SISO

% reduced system CMS_1
M_red_CMS_1 = T_CMS_1'*M*T_CMS_1;
K_red_CMS_1 = T_CMS_1'*K*T_CMS_1;
D_damp_red_CMS_1 = T_CMS_1'*D_damp*T_CMS_1;

ndof_red_CMS_1 = rank(M_red_CMS_1);

A_red_CMS_1 = [zeros(ndof_red_CMS_1) eye(ndof_red_CMS_1);-M_red_CMS_1\...
    K_red_CMS_1 -M_red_CMS_1\D_damp_red_CMS_1]; %M\K equals inv(M)*K
B_red_CMS_1 = [zeros(ndof_red_CMS_1);inv(M_red_CMS_1)];
C_red_CMS_1 = [eye(ndof_red_CMS_1) zeros(ndof_red_CMS_1)];
D_red_CMS_1 = zeros(ndof_red_CMS_1);

sys_red_CMS_1 = ss(A_red_CMS_1,B_red_CMS_1,C_red_CMS_1,...
    D_red_CMS_1); %create MIMO state space system
sysss_red_CMS_1 = sss(sys_red_CMS_1); %convert into sparse state space system
sysss_red_CMS_1_trunc = truncate(sysss_red_CMS_1,...
    output_DOF,input_DOF); %reduce to SISO

% reduced system CMS_2
M_red_CMS_2 = T_CMS_2'*M*T_CMS_2;
K_red_CMS_2 = T_CMS_2'*K*T_CMS_2;
D_damp_red_CMS_2 = T_CMS_2'*D_damp*T_CMS_2;

ndof_red_CMS_2 = rank(M_red_CMS_2);

A_red_CMS_2 = [zeros(ndof_red_CMS_2) eye(ndof_red_CMS_2);-M_red_CMS_2\...
    K_red_CMS_2 -M_red_CMS_2\D_damp_red_CMS_2]; %M\K equals inv(M)*K
B_red_CMS_2 = [zeros(ndof_red_CMS_2);inv(M_red_CMS_2)];
C_red_CMS_2 = [eye(ndof_red_CMS_2) zeros(ndof_red_CMS_2)];
D_red_CMS_2 = zeros(ndof_red_CMS_2);

sys_red_CMS_2 = ss(A_red_CMS_2,B_red_CMS_2,C_red_CMS_2,...
    D_red_CMS_2); %create MIMO state space system
sysss_red_CMS_2 = sss(sys_red_CMS_2); %convert into sparse state space system
sysss_red_CMS_2_trunc = truncate(sysss_red_CMS_2,...
    output_DOF,input_DOF); %reduce to SISO

%% calculate values for FRF

% magnitude and phase of bode plot
[mag_full, phase_full] = bode(sysss_full_trunc,sampling_frequ);
[mag_red_dyn, phase_red_dyn] = bode(sysss_red_dyn_trunc,sampling_frequ);
[mag_red_CMS_1, phase_red_CMS_1] = bode(sysss_red_CMS_1_trunc,sampling_frequ);
[mag_red_CMS_2, phase_red_CMS_2] = bode(sysss_red_CMS_2_trunc,sampling_frequ);

% relative error in magnitude
mag_rel_error_red_dyn = abs(mag_red_dyn-mag_full)./mag_full;
mag_rel_error_red_CMS_1 = abs(mag_red_CMS_1-mag_full)./mag_full;
mag_rel_error_red_CMS_2 = abs(mag_red_CMS_2-mag_full)./mag_full;

%% plot FRFs and difference

figure;
loglog(sampling_frequ,squeeze(mag_full),'k-')
hold on
loglog(sampling_frequ,squeeze(mag_red_dyn),'g-*')
loglog(sampling_frequ,squeeze(mag_red_CMS_1),'b--')
loglog(sampling_frequ,squeeze(mag_red_CMS_2),'r-.')
dynamic_str = strcat('dynamic w=',num2str(frequ_dynamic_reduction));
legend('original',{'full',dynamic_str,'CB 20 modes','CB 40 modes'},...
    'Location','NorthEast')
title('FRF comparison')
xlabel('Frequency [rad]')
ylabel('Magnitude of FRF')

figure;
loglog(sampling_frequ,squeeze(mag_rel_error_red_dyn),'g-*')
hold on
loglog(sampling_frequ,squeeze(mag_rel_error_red_CMS_1),'b--')
loglog(sampling_frequ,squeeze(mag_rel_error_red_CMS_2),'r-.')
legend('original',{dynamic_str,'CB 20 modes','CB 40 modes'},...
    'Location','NorthEast')
title('Rel. error comparison of FRF')
xlabel('Frequency [rad]')
ylabel('Rel. error [-]')

% h = figure;
% bode(sysss_full_trunc,'k-');
% hold on
% bode(sysss_red_dyn_trunc,'g-*');
% bode(sysss_red_CMS_1_trunc,'b--');
% bode(sysss_red_CMS_2_trunc,'r-.');
% legend('original',{'full','dynamic','CB_20_modes','CB_40_modes'},...
%     'Location','NorthEast')