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

% apply boundary conditions
% BC_set = [1:12];
% 
% K = K;
% K(BC_set,:) = [];
% K(:,BC_set) = [];
% M = M;
% M(BC_set,:) = [];
% M(:,BC_set) = [];
% D = D_damp;
% D(BC_set,:) = [];
% D(:,BC_set) = [];

% main parameters for methods
master_set = [1:12]; %for CMS and Dynamic reduction
modes_to_keep_CMS_1 = 1:20;

% parameters for FRF comparison
output_DOF = 1;
input_DOF = 1;
sampling_frequ = 7e3:.1:5e4; %in rad/s

%% FRF comparison

% full system
sys_full = exp_ss_model( M, D_damp, K );

sysss_full = sss(sys_full); %convert into sparse state space system
sysss_full_trunc = truncate(sysss_full,output_DOF,input_DOF); %reduce to SISO

% reduced system CMS_1
[sys_red_CMS_1, M_red_CMS_1, D_damp_red_CMS_1, ...
    K_red_CMS_1, T_CMS_1] = CMS_reduction(...
    'CB', modes_to_keep_CMS_1, M, D_damp, K, master_set);

sysss_red_CMS_1 = sss(sys_red_CMS_1); %convert into sparse state space system
sysss_red_CMS_1_trunc = truncate(sysss_red_CMS_1,...
    output_DOF,input_DOF); %reduce to SISO

%% calculate values for FRF

% magnitude and phase of bode plot
[mag_full, phase_full, frequ_rad] = bode(sysss_full_trunc,sampling_frequ);
[mag_red_CMS_1, phase_red_CMS_1] = bode(sysss_red_CMS_1_trunc,sampling_frequ);

% relative error in magnitude
mag_rel_error_red_CMS_1 = abs(mag_red_CMS_1-mag_full)./mag_full;

%% plot FRFs and difference

frequ_hz = frequ_rad/(2*pi);

figure;
loglog(frequ_hz,squeeze(mag_full),'k-')
hold on
loglog(frequ_hz,squeeze(mag_red_CMS_1),'b--')
legend('original',{'full','CB'},'Location','NorthEast')
title('FRF comparison')
xlabel('Frequency [Hz]')
ylabel('Magnitude of FRF')

figure;
loglog(frequ_hz,squeeze(mag_rel_error_red_CMS_1),'b--')
legend('original',{'CB',},'Location','NorthEast')
title('Rel. error comparison of FRF')
xlabel('Frequency [Hz]')
ylabel('Rel. error [-]')



