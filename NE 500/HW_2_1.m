clear all, close all, clc
%% #################### HW-2.1 ##################################
% Units: P = KPa, T = K 
% h = J/kg, s = J/kg-k, v = m^3/kg

P_1 = 6880 ; %(kPa) 68.8 bar
h_1 = 2774.12;

% at point 2
h_2f = 811.58;
h_2g = 2785.98;
h_2fg = h_2g - h_2f; 
s_1 = 5.82;
s_2f = 2.24;
s_2g = 6.49;
s_2fg = s_2g - s_2f;

eff_T = 0.9; % Turbine efficiency

x_2s = (s_1 - s_2f) /s_2fg;
h_2s = h_2f + x_2s *h_2fg;
W_ts_hp = h_1 - h_2s;
W_ta_hp = W_ts_hp *eff_T;
h_2a = h_1-W_ta_hp;
x_2a = (h_2a - h_2f)/h_2fg;
%m2=(1-x_2a);
mR = 1 - x_2s;

% at point 3
P_2 = 3;%% (kPa) 0.03 (bar)
h_3 = 2785.98; %saturated

s_3 = 6.49;
s_5s = s_3;

% at point 5%%
s_5f = 0.35;
s_5g = 8.57;
s_5fg = s_5g - s_5f;
x_5s = (s_5s - s_5f)/s_5fg;

h_5f = 100.99;
h_5g = 2544.99;
h_5fg = h_5g - h_5f;
h_5s = h_5f + x_5s *h_5fg;

% x_5a = (h_5s - h_5f) /h_5fg;
W_ts_lp_m3 = (h_3-h_5s); % m3 = m1 - m2
W_ts_lp_m31 = W_ts_lp_m3 *x_2s;
W_ts_lp = x_2a *W_ts_lp_m3; %% m3 = m1 - m2
W_ta_lp = eff_T *W_ts_lp;

% at point 6
v_6 = 0.001;
P_7 = 1280;
P_6 = 3;
eff_P = 0.8;
W_cps = v_6*(P_7-P_6)*x_2a;
W_cps1 = v_6*(P_7-P_6)*x_2s;
W_cpa = W_cps/eff_P;

% at point 8
v_8 = 0.0011;
P_9 = 6880;
P_8 = 1280;
W_fps = v_8 *(P_9 - P_8);
W_fpa = W_fps/eff_P;
h_6 = h_5f;
h_7a = h_6 + W_fpa;
h_7s = h_6 + W_fps;
h_4s = 811.57;
h_8a = h_4s *(1-x_2a) + h_7a *x_2a;
h_8s = h_4s *(1-x_2s) + h_7s *x_2s;
h_9s = h_8s + W_fps;
h_9a = h_8a + W_fpa;
Q_h = h_1 - h_9a;
Q_hs = h_1 - h_9s;

%% Efficiency
% Actual Efficiency
eff_a =(W_ta_hp + W_ta_lp - W_fpa - W_cpa)/Q_h;

% Ideal Efficiency
eff_s =(W_ts_hp + W_ts_lp_m31- W_fps - W_cps1)/Q_hs;


%% Lost Work
T_env = 300;  % 27 dC

%### TBC

