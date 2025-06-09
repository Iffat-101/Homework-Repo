clear all, close all, clc
%% #################### HW-2.5 ##################################
%% calculation, section a
    %% No. of available cores
N_c = 8730112;  % total no. of cores
a_v = 0.25;     % availability factor
N_a = N_c *a_v; % cores available

    %% No. of processable nodes
n_o = 2048;     % processable nodes /core
N_d = N_a *n_o;  % processable nodes, total, volume

    %% Total Volume of bubble
v_f = 0.02;     % const. factor of bubble
V_bt = N_d *v_f; % total bubble volume

    %% Number of bubbles
n_d = 10;       % mesh resln @ bubble dia
% n_bb = V_bt *6 /(pi *n_d^3);
syms vr         % declaring variable
exp = V_bt == vr *pi*n_d^3 /6;
n_bb = double(solve(exp, vr));

