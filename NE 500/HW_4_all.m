clear all, close all, clc
%% #################### HW_4 ##################################
% Units: P = Psia, T = 0F 
% h = BTU/lb, s = BTU/lb-F, v = ft^3/lb

D_CF = 1.5/12; D_CB = 7.8/12;

A_CF = pi * D_CF^2/4;
A_CB = pi * D_CB^2/4;
A_ds = A_CB /A_CF;
z_1 = 0.5/12; t_CB = 1.3/12;
z_2 = z_1 + t_CB;

T_z3 = -319; % Ambient temperature
T_z1 = 1981; % Copper's MP
k_cu = 223; k_W = 94;

Q_dt = 310 *3413;
q_ddt_0 = Q_dt /A_CB; % heat flux at beginning


T_z2 = T_z1 - q_ddt_0 /k_cu *(z_2 - z_1);

dl_z32 = (T_z2 - T_z3) * k_cu / q_ddt_0 * 1/A_ds;
dl_z32_in = dl_z32 *12; 

% For Tungsten as limitting factor
T_z0 = 6170;
T_z1_W =  T_z0 - q_ddt_0 /k_W *z_1;
T_z2_W = T_z1_W - q_ddt_0 /k_cu *(z_2 - z_1);
dl_z32_W = (T_z2_W - T_z3) * k_cu / q_ddt_0 * 1/A_ds;
dl_z32_W_in = dl_z32_W *12;