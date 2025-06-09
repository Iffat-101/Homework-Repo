clear all, close all, clc
%% #################### HW-1.3 ##################################
%% Givens and Table Readings
P_1 = 1250; 

P_1d = P_1; T_1d = 572.46;
h_1d = 578.76; s_1d = 0.78;

P_2d = P_1; T_2d = T_1d; 
h_2d = 1181.89; s_2d = 1.36;

P_2 = P_2d; T_2 = 690;
h_2 = 1299.09; s_2 = 1.47;

P_4d = 430; T_4d = 451.77;
h_4d = 1205.17;  s_4d = 1.48;

P_4 = P_4d; T_4 = T_4d;
h_4 = 432.19; s_4 = 0.63;   v_4 = 0.0195;
        % this is identical with state 9

P_3 = 250; T_3 = 400.98;
        s_3 = s_2;  % isentropic expansion
 
P_3g = P_3; T_3g = 400.98;
h_3g = 1201.59; s_3g = 1.53;

P_3f = P_3; T_3f = T_3;
h_3f = 376.16; s_3f = 0.57;

P_5 = 250; T_5 = 451.77; % Big assumption
h_5 = 1235.06; s_5 = 1.56; % Superheated

P_6 = 2.5; T_6 = 134.38;
        s_6 = s_5; % isentropic expansion

P_6g = P_6; T_6g = T_6;
h_6g = 1119.26; s_6g = 1.9;

P_7 = P_6; T_7 = T_6;
h_7 = 102.36; s_7 = 0.189; v_7 = 0.016;

P_8 = 430;

P_9 = P_4; T_9 = T_4;
h_9 = h_4; s_9 = s_4; v_9 = v_4;

Cf = 5.4;

%% Calculation
% Steam quality at 3
x_3 = (s_3 - s_3f) /(s_3g - s_3f);

% Enthalpy at 3
h_3 = h_3f + x_3 *(h_3g - h_3f);

% Steam Quality at 6
x_6 = (s_6 - s_7) /(s_6g - s_7);

% Enthalpy at 6
h_6 = h_7 + x_6 *(h_6g - h_7);

% Heat Transfer in Reheater
    % this is where mass ratio mR is found
mR = (h_5 - h_3) /((h_2 + h_5) - (h_4 + h_3));

% HP Turbine, Output 
W_t_hp =  h_2 - h_3 *(1 - mR);

% LP Turbine, Output
W_t_lp =  (h_5 - h_6) *(1 - mR);

% Condenser Cooling
W_l = (h_6 - h_7) * (1 - mR);


% Condenser Pump, input
    % enthalpy at 8
h_8 = h_7 + v_7 *(P_8 - P_7) /Cf;
W_cp = (1 - mR) * (h_8 - h_7);

% Feed pump, input
    % enthalpy at 1
h_1 = h_9 + v_9 * (P_1 - P_9) /Cf;
W_fp = h_1 - h_9;

% Heating in Boiler
Q_h = h_2 - h_1;

% Efficiency
W_ot = W_t_hp + W_t_lp;
W_in = W_fp + W_cp;
W_nt = W_ot - W_in;

eff = W_nt /Q_h;





