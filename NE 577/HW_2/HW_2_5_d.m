clear all, close all, clc
%% #################### HW-2.5-1d ##################################
    %% No. of processable nodes
Cp = 3.25 *10^12;       % computational power
Cp_Nd = Cp /10^6;       % comptnl power per node
Cp_mx = 1.102 *10^18;   % Max compt. power
N_d =  Cp_mx /Cp_Nd;    % N_d for prob. d

    %% Total volume of bubble
v_f = 0.02;     % constitute factor of bubble
V_bt = N_d *v_f; % total bubble volume

 %% Number of bubbles
n_d = 10;       % mesh resln @ bubble dia
% n_bb = V_bt *6 /(pi *n_d^3);
syms vr         % declaring variable
exp = V_bt == vr *pi*n_d^3 /6;
n_bb = double(solve(exp, vr));
