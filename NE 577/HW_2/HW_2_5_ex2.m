clear all, close all, clc
%% #################### HW-2.5-1d ##################################
    %% No. of processable nodes
% for fastest computer at June 2013
Cp = 3.25 *10^12;
Cp_Nd = Cp /10^6;
Cp_mx = 33.862 *10^15;   % Top500, @ '13
N_d =  Cp_mx /Cp_Nd;    % N_d for prob. d

    %% Total volume of bubble
v_f = 0.02;     % constitute factor of bubble
V_bt = N_d *v_f; % total bubble volume

    %% Number of bubbles
n_d = 20;       % mesh resln @ bubble dia
% n_bb = V_bt *6 /(pi *n_d^3);
syms vr         % declaring variable
exp = V_bt == vr *pi *n_d^3 /6;
n_bb = double(solve(exp, vr));

%% Part 2 ####################
CP_mx_13 = 33.862 *10^15;
CP_mx_23 = 1.102 *10^18;
gr = CP_mx_23 /CP_mx_13; % growth

N_c_13 = 3120000;
N_c_23 = 8730112;
gr_cr =  N_c_23 /N_c_13;  % core growth
N_c_33 = gr_cr * N_c_23;

    %% So, No. of processable nodes
Cp = 3.25 *10^12;
Cp_Nd = Cp /10^6;
CP_mx_33 = gr * CP_mx_23;
N_d2 = CP_mx_33 /Cp_Nd;

    %% Total volume of bubble
v_f2 = 0.02;     % constitute factor of bubble
V_bt2 = N_d2 *v_f2; % total bubble volume

    %% Number of bubbles
n_d2 = 20;       % mesh resln @ bubble dia
syms vr2         % declaring variable
exp = V_bt2 == vr2 *pi *n_d2^3 /6;
n_bb2 = double(solve(exp, vr2));


