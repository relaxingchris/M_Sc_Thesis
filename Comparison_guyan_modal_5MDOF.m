close all
clear all
clc

% Comparison of GUYAN reduction, modal reduction and
% Component-Mode-Synthesis with an example beam

addpath('/Users/LaKrizz/Documents/MATLAB/Master_Thesis/model_data')

% read mass and stiffness matrix from file
%M = generate_matrix_from_mtx_file('Bar_more_detailed_v2_MASS2.mtx');
M = generate_matrix_from_mtx_file('Beam_3_element_rectangle_COORDINATE_MASS2.mtx');

%K = generate_matrix_from_mtx_file('Bar_more_detailed_v2_STIF2.mtx');
K = generate_matrix_from_mtx_file('Beam_3_element_rectangle_COORDINATE_STIF2.mtx');


% generate Damping matrix
alpha = 0.2;
beta = 0;
D_damp = alpha*M + beta*K;

%DOF_interest = [546];
DOF_interest = [15*3];

% generate Input vector
Inp_vec = zeros(rank(M),1);
%DOF_input = [546];
DOF_input = [15*3];
for i=DOF_input
    Inp_vec(i) = 1;
end

%% choose which reduction methods to use
guyan_red = 'no';
master_DOF_guy = [1:(9*3) (90*3+1):(99*3) (180*3+1):(189*3)];

modal_red = 'yes';
modes_to_keep = 1:20;

%% solve full system
[eig_vec_full, eig_val_full] = eig_vec_and_eig_val_hz(M, K);

% COMMENTED FOR DEBUGGING!

if 1
    
    % set up state space system
    %A_full = [zeros(rank(M)) eye(rank(M));-inv(M)*K zeros(rank(M))];
    % with damping:
    A_full = [zeros(rank(M)) eye(rank(M));-inv(M)*K -inv(M)*D_damp];
    B_full = [zeros(rank(M));inv(M)];
    B_full = B_full*Inp_vec; %reduce system to only have the defined inputs
    C_full = [zeros(rank(M)) zeros(rank(M))];
    C_full(DOF_interest,DOF_interest) = 1; %reduce output to only DOF_interest
    D_full = zeros(rank(M));
    D_full = D_full*Inp_vec; %reduce system to only have the defined inputs

    sys_full = ss(A_full,B_full,C_full,D_full);

    eig_val = eig(sys_full);
    eig_val = sort(eig_val);
    eig_val_hz = eig_val/(2*pi);
    eig_val_hz = abs(eig_val_hz);

    min_freq = 1; %Hz
    max_frequ = 1000; %Hz
    frequ_rad = min_freq*2*pi:1:max_frequ*2*pi; %Frequency for bode plot in [rad/s]

    % calculate magnitude and phase of frequency response function
    [mag, phase, frequ] =  bode(sys_full,frequ_rad);

    % change from [rad/s] to [Hz]
    frequ_hz = frequ/(2*pi);

    % select FRF of DOF where load is applied
    mag_point_of_int_full = squeeze(mag(DOF_interest(1),1,:));
    phase_point_of_int_full = squeeze(phase(DOF_interest(1),1,:));

    figure
    subplot(2,1,1)
    plot(frequ_hz,20*log10(mag_point_of_int_full))
    xlabel('Frequency [Hz]')
    ylabel('Magnitue [dB]')
    subplot(2,1,2)
    plot(frequ_hz,phase_point_of_int_full)
    xlabel('Frequency [Hz]')
    ylabel('Phase [deg]')
end


%% static (GUYAN) reduction

if (strcmpi(guyan_red,'yes'))
    % set up slave DOF
    slave_DOF = 1:rank(M);
    slave_DOF(ismember(slave_DOF,master_DOF_guy)) = [];
    
    % sort matrices
    M_mm = M(master_DOF_guy,master_DOF_guy); 
    M_ms = M(master_DOF_guy,slave_DOF);
    M_sm = M(slave_DOF,master_DOF_guy);
    M_ss = M(slave_DOF,slave_DOF);

    M_sort = [M_mm M_ms;
        M_sm M_ss];

    K_mm = K(master_DOF_guy,master_DOF_guy); 
    K_ms = K(master_DOF_guy,slave_DOF);
    K_sm = K(slave_DOF,master_DOF_guy);
    K_ss = K(slave_DOF,slave_DOF);

    K_sort = [K_mm K_ms;
        K_sm K_ss];
    
    %create transformation matrix
    T_guy = [eye(size(master_DOF_guy,2));
        -inv(K_ss)*K_sm];

    %transform system
    [M_red_guy, K_red_gyu] = transform_system(M_sort, K_sort, T_guy);
    
    %solve reduced system
    [eig_vec_red_guy, eig_val_red_guy] = eig_vec_and_eig_val_hz(...
    M_red_guy, K_red_gyu);

    % compare eigenfrequencies
    NFD_guy_full = nat_frequ_diff_criterion( eig_val_full,...
    eig_val_red_guy, master_DOF_guy);
    
    % modMAC
    eig_vec_red_guy_retrans = T_guy*eig_vec_red_guy;

    mac_guy_full = MAC_mass_mod(eig_vec_full, eig_vec_red_guy_retrans, M,...
    master_DOF_guy, 'real');
end

%% Modal reduction

if (strcmpi(modal_red,'yes'))
    %create transformation matrix
    T_mod = eig_vec_full(:,modes_to_keep);

    %transform system
    [M_red_mod, K_red_mod] = transform_system(M, K, T_mod);
    
    %solve reduced system
    [eig_vec_red_mod, eig_val_red_mod] = eig_vec_and_eig_val_hz(...
    M_red_mod, K_red_mod);
    
    %compare eigenfrequencies
    NFD_mod_full = nat_frequ_diff_criterion(eig_val_full,...
    eig_val_red_mod, modes_to_keep);
    
    %make eigenvectors orthogonal
    %eig_vec_red_mod = orth(eig_vec_red_mod);
    
    %modMAC
    eig_vec_red_mod_retrans = T_mod*eig_vec_red_mod;
    
    mac_mod_full = MAC_mass_mod(eig_vec_full, eig_vec_red_mod_retrans, M,...
    modes_to_keep, 'real');
end




