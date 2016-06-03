% clear all
% close all
% clc

% read mass and stiffness matrix from file
%M = generate_matrix_from_mtx_file('Beam_3_element_FORMAT_COORDINATE_MASS2.mtx');
%K = generate_matrix_from_mtx_file('Beam_3_element_FORMAT_COORDINATE_STIF2.mtx');

% rectangular cross section
M = generate_matrix_from_mtx_file('Beam_3_element_rectangle_COORDINATE_MASS2.mtx');
K = generate_matrix_from_mtx_file('Beam_3_element_rectangle_COORDINATE_STIF2.mtx');

% remove rigid body modes
%K = K([13:48],[13:48]);
%M = M([13:48],[13:48]);

% solve full system
[eig_vec_full, eig_val_full] = eig_vec_and_eig_val_hz(M, K);

% choose which reduction methods to use
guyan_red = 'yes';
master_DOF_guy = [1:12 37:48];

modal_red = 'no';
modes_to_keep = 1:20;

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







