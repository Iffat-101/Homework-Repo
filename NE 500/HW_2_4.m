clear all, close all, clc
%% #################### HW-2.4 ##################################
%% Givens and Table Readings
% Units: P = Psi, T = 0F 
% h = BTU/lb, s = BTU/lb-F v = lb/ft^3

P_1 = 780; T_1 = nan;

P_2 = 780; T_2 = 800;
h_2 = 1399.56; s_2 = 1.6;

P_2g = P_2; T_2g = 515.36;

P_4 = 3; T_4 = 141.42;
         s_4 = s_2;

% 5 = 4f 
P_5 = P_4; T_5 = T_4;
h_5 = 109.39; s_5 = 0.201; v_5 = 0.0163;

P_4g = P_4; T_4g = T_4;
h_4g = 1122.19; s_4g = 1.886;

cf = 5.4; % BTU/ft^3.psi

%% Table Readings, from Calculated
% Temperature midway between boiler-saturation and condenser
T_st_bl = T_2g; T_cd = T_4;
T_md = (T_st_bl + T_cd) /2;

% Now, data will be read for T_3 = T_md = 328.39;
T_3 = 328.39; P_3 = 100.8;
              s_3 = s_2; % isentropic expansion

P_3f = P_3; T_3f = T_3;
h_3f = 299.18; s_3f = 0.475; v_3f = 0.0177;

P_3g = P_3; T_3g = T_3;
h_3g = 1187.63; s_3g = 1.603;

P_6 = P_3;

P_7 = P_3; T_7 = 328.39;
h_7 = h_3f; s_7 = s_3f; v_7 = v_3f;

%% Calculation
% High pressure tap
    % steam quality at 3
    s_3 = s_2; 
    x_3 = (s_3 - s_3f) /(s_3g - s_3f);

    % enthalpy at 3
    h_3 = h_3f + x_3 * (h_3g - h_3f);

% Turbine outlet
    % quality at 4
    s_4 = s_2; % isentropic expansion
    x_4 = (s_4 - s_5) /(s_4g - s_5);

    % enthalpy at 4
    h_4 = h_5 + x_4 *(h_4g - h_5);

% Condenser Pump Work, Input
W_cp_m3 = v_5 *(P_6 - P_5) /cf;
    % this is normalized by m3 = m1 - m2
    
    % enthalpy at 6
    h_6 = h_5 + W_cp_m3;

    % Feedwater Heater
    % by first law, we get mR = m1/m2
    mR = (h_7 - h_6) /(h_3 - h_6); %# check it 
    
    % Condenser pump work
    W_cp_m1 = (1 - mR) * W_cp_m3;
        % this is normalized by m1

% Feed Pump Work, Input
W_fp_m1 = v_7 *(P_1 - P_7) /cf;

    % enthalpy at 1
    h_1 = h_7 + W_fp_m1;

% Boiler heat, input
Q_h_m1 = h_2 - h_1;

% Turbine Work, Output
W_t_m1 = h_2 - h_3 * mR - h_4 *(1-mR);


% Cycle Efficiency
    % Net work
    W_nt = W_t_m1 - (W_cp_m1 + W_fp_m1);

    % efficiency;
    eff = W_nt / Q_h_m1 *100;
 



