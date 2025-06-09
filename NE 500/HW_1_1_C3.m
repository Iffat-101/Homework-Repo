clear all, close all, clc
%% #################### HW-1.1.C3 ##################################
%% Givens and Table Readings
% Units: P = Psi, T = 0F 
% h = BTU/lb, s = BTU/lb-F v = ft^3/lb

T_5 = 502.5; P_5 = 696.18;
h_5 = 1201.99; s_5 = 1.43; v_5 = 0.66; x_5 = 1;

T_4 = 502.5; P_4 = 696.18;
h_4 = 490.89;   s_4 = 0.69;  v_4 = 0.02;   x_4 = 0;

T_1 = 705; P_1 = 696.18;  % P_1 = 4.8 MPa 
h_1 = 1348.05; s_1 = 1.57; v_1 = 0.92; x_1 = 1;

T_2 = 90; P_2 = 0.699;
s_2 = s_1; % isentropic expansion

T_2g = 90; P_2g = 0.699;
h_2g = 1100.43; s_2g = 2.0;  v_2g = 467.5; x_2g = 1;

T_3 = 90; P_3 = 0.699;
h_3 = 58.05;   s_3 = 0.112;    v_3 = 0.016;   x_3 = 0;

P_3d = P_4;

Cf = 5.4;   % for British unit 
Cf2 = 2.326; % 1 BTU/lb -> KJ/kg

%% Calculation
% steam quality at 2
s_2 = s_1; % isentropic expansion
x_2 = (s_2 - s_3) /(s_2g - s_3);

% enthalpy at 2
h_2 = h_3 + x_2 *(h_2g - h_3);

% Turbine work
W_t = h_1 - h_2;

% Pump work
W_p = v_3 *(P_3d - P_3) /Cf; % this is -ve, as given to system

% Enthalpy at 3'
h_3d = h_3 + W_p;

% Heating in Boiler
Q_h = h_1 - h_3d;

% Cooling in Condenser
Q_L = h_2 - h_3;

% Efficiency
eff = (W_t - W_p) /Q_h;

% Steam Rate
W_nt = (W_t - W_p) * Cf2; % KJ/Kg 
s_rt = 3600 /(W_nt);


%% Heat Added in legs
% Heat added in leg 3-4
h_3d_4 = h_4 - h_3d;

% in leg 4-1
h_4_1 = h_1 - h_4;
