clear all, close all, clc
%% #################### HW-1.1.C2 ##################################
%% Givens and Table Readings
% Units: P = Psi, T = 0F 
% h = BTU/lb, s = BTU/lb-F v = ft^3/lb

T_1 = 695; P_1 = 2991.9; 
h_1 = 1018.31; s_1 = 1.16;    v_1 = 0.085;    x_1 = 1;

T_2g = 90; P_2g = 0.699;
h_2g = 1100.43; s_2g = 2.0;                 x_2g = 1;

T_3f = 90; P_3f = 0.699;
h_3f = 58.05;   s_3f = 0.112;    v_3f = 0.016;   x_3f = 0;

T_4 = 695; P_4 = P_1;
h_4 = 801.35;   s_4 = 0.97;    v_4 = 0.034;   x_4 = 0;

T_2 = 90; P_2 = P_2g;
s_2 = s_1;  % isentropic process

Cf = 5.4;
Cf2 = 2.326; % 1 BTU/lb -> KJ/kg

%% Calculation
% steam quality at 2
s_2 = s_1; % isentropic expansion
x_2 = (s_2 - s_3f) /(s_2g - s_3f);

% enthalpy at 2
h_2 = h_3f + x_2 *(h_2g - h_3f);

% Quality at 3
s_3 = s_4;
x_3 = (s_3 - s_3f) /(s_2g - s_3f);

% Enthalpy at 3
h_3 = h_3f + x_3 *(h_2g - h_3f);

% Turbine work, Output
W_t = h_1 - h_2;

% Pump work, Input
W_p = h_4 - h_3;

% Heating in Boiler
Q_h = h_1 - h_4;

% Cooling in Condenser
Q_L = h_2 - h_3;

% Efficiency
eff = (W_t - W_p) /Q_h;

% Steam Rate
W_nt = (W_t - W_p) * Cf2; % KJ/Kg 
s_rt = 3600 /(W_nt);

%% Heat Added in legs
% Heat added in leg 3-4
h_3d_4 = h_4 - h_3;

% in leg 4-1
h_4_1 = h_1 - h_4;

