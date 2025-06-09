clear all, close all, clc
%% #################### HW-2.3 ##################################
%% Givens and Table Readings
% Units: P = Psi, T = 0F 
% h = BTU/lb, s = BTU/lb-F, v = ft^3/lb

P_2 = 1220;   T_2 = 569.36;
h_2 = 1152.86; s_2 = 1.336;

s_3f = 0.72;
s_3g = 1.41;
s_3fg = s_3g - s_3f;

x_3 = (s_2 - s_3f)/s_3fg;
h_3f = 516.73;
h_3g = 1198.1;
h_3fg = h_3g - h_3f;
h_3 = h_3f + x_3 *h_3fg;

% at point 4
P_4 = 350;
s_4 = 1.37;
s_4f = 0.61; s_4g = 1.497;
s_4fg = s_4g-s_4f;
x_4 = (s_4-s_4f)/s_4fg;
m5_m4 = (1 - x_4);
h_4f = 409.84;
h_4g = 1204.47;
h_4fg = h_4g - h_4f;
h_4 = h_4f + x_4 *h_4fg ;

% at point 7
h_7 = 1196;

% point at 8 
s_7 = 1.4;

s_8f = 0.2;
s_8g = 1.88;
s_8fg = s_8g - s_8f ;
x_8 = (s_7-s_8f) /s_8fg ;
h_8f = 109.39;
h_8g = 1122.19;
h_8fg = h_8g - h_8f;
h_8 = h_8f + x_8 *h_8fg;

% at point 6%
h_6 = 576.62;

% at point 1 %
h_1 = 369 ;
% h_2 = 1183.29 ;

% at point 14
s_14f = 0.43;
s_14g = 1.64;
s_14fg = s_14g - s_14f;
x_14 = (s_7 - s_14f) /s_14fg;
h_14f = 267.66;
h_14g = 1179.36;
h_14fg = h_14g - h_14f;

h_14= h_14f + x_14 * h_14fg;

% At point 10
h_9 = 109.39  ;
v_9 = 0.0162;
P_10 = 65;
P_9 = 3;
CF = 5.4;
W_cp = v_9 *(P_10 - P_9)/CF;
h_10 = h_9 + W_cp;

% Enthalpy at point 5
h_5 = 409.85; 

% at point 11
h_11 = 267.66;  v_11 = 0.017;
P_12 = 1220;
P_11 = 65;
W_fp = v_11*(P_12 - P_11) /CF;
h_12 = h_11 + W_fp;

% at point 13
h_13 = 516.72; s_13 = 0.718;

% mass ratios
x_4m = 1- x_4; % moisture content at 4
mr3 = (h_1 - h_12)/(h_3 - h_13);

% %     syms vr4
% exp = vr4/x_4m *h_4 + (1 - vr4 /x_4m - mr3) *h_2 == vr4 *h_5 + (1- vr4/x_4m - mr3) *h_6 + ...
% vr4 *(1 /x_4m - 1) *h_7;
% mr4 = double(solve(exp, vr4));  % shows into fraction

mr2 = 0.127 *(1 - mr3);
mr4 = x_4m *(1- mr2 - mr3);
mr5 = ( h_11 - (1 -mr2 - mr3 -mr4) *h_10 - mr2 *h_6 - mr4 *h_5 - mr3 *h_13 ) /(h_14 - h_10);

% Boiler heat input
Q_h_m1 = h_2 - h_1;

% High pressure turbine work
mr3_ds = 1 - mr2 - mr3;
mr2_ds = 1 -mr2;
W_t_hp_m1 = h_2 * mr2_ds - h_3 *mr3 - mr3_ds *h_5;

% Low pressure turbine work
    mr4_ds = 1- mr2 - mr3 - mr4; 
    mr5_ds = 1- mr2 - mr3 - mr4 - mr5;
    W_t_lp_m1 = mr4_ds * h_7 - mr5 *h_14 - mr5_ds *h_8;

% Condenser Pump Work
W_cp_m1 = mr5_ds *W_cp;
 
% Feed Pump Work
W_fp_m1 = W_fp;

% efficiency
eff = (W_t_lp_m1 + W_t_lp_m1 - W_cp_m1 - W_fp_m1) / Q_h_m1 * 100; 

%% Test Value
val = x_4 *(h_7 - h_4g) /(h_2 - h_6);
val2 = x_4 *(h_7 - h_4) /(h_2 - h_6);
