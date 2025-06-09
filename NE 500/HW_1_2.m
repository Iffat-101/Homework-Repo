clear all, close all, clc
%% #################### HW-1.2 ##################################
%% Givens and Table Readings
% Units: P = Psi, T = 0F 
% h = BTU/lb, s = BTU/lb-F v = ft^3/lb
P_1 = 1250;

P_1d = 1250; T_1d = 572.46;
h_1d = 578.76; s_1d = 0.778;    v_1d = 0.023;  x_1d = 0;

T_2d = T_1d; P_2d = P_1d;
h_2d = 1181.89; s_2d = 1.36;  v_2d = 0.345;  x_2d = 1;

T_2 = 690; P_2 = P_2d;
h_2 = 1299.09; s_2 = 1.47;  v_2 = 0.46; x_2 = 1;

% T_3 , P_3 = unknown
T_4 = 690;

P_5 = 2.5; T_5 = 134.38;

P_5g = P_5; T_5g = T_5;
h_5g = 1119.26; s_5g = 1.9; v_5g = 140.9;  x_5g = 1;

P_6 = P_5; T_6 = T_5;
h_6 = 102.37; s_6 = 0.189; v_6 = 0.016; x_6 = 0;

Cf = 5.4;

%% Design Considerations
P_3 =   [150,     200,     250,     300,     350,     400 ];
T_3 =   [358.43,  381.81,  400.98,  417.36,  431.75,  444.62 ];
h_3g =  [1194.49, 1198.80, 1201.58, 1203.38, 1204.47, 1205.04 ]; % Sat vapor
s_3g =  [1.57,    1.54,    1.52,    1.51,    1.49,    1.48 ];       
h_3f =  [330.68,  355.53,  376.16,  393.99,  409.84,  424.17];   % Sat liq
s_3f =  [0.514,   0.543,   0.568,   0.58,    0.61,    0.62];

P_4  =  [150,     200,     250,     300,     350,     400 ];
h_4 =   [1371.69, 1368.91, 1366.09, 1363.24  1360.34  1357.39];
s_4 =   [1.75,    1.71,    1.69,    1.67,    1.65,    1.63];

%% Table Readings from calculated findings
T_4 = 690; 

%% Calculation

for i = 1:1:6
    % steam quality at 3
    s_3(i) = s_2; % isentropic expansion
    x_3(i) = (s_3(i) - s_3f(i)) /(s_3g(i) - s_3f(i));

    % Enthalpy at 3
    h_3(i) = h_3f(i) + x_3(i) * (h_3g(i) - h_3f(i));
    
    % HP Turbine Work
    W_t_hp(i) = h_2 - h_3(i);
    
    % Steam Quality at point 5
    s_5(i) = s_4(i); % isentropic expansion
    x_5(i) = (s_5(i) - s_6) /(s_5g - s_6);

    % Enthalpy at 5
    h_5(i) = h_6 + x_5(i) * (h_5g - h_6);
    
    % LP Turbine Work
    W_t_lp(i) = h_4(i) - h_5(i);
    
    % Pump work
    W_p = v_6 * (P_1 - P_6) /Cf; % this is -ve, as given to system
    
    
    % Enthalpy at 1  
    h_1 = W_p + h_6;
    
    % Heating in Boiler
    Q_h(i) = h_4(i) + h_2 - (h_1 + h_3(i));
    
    % Cooling in Condenser
    Q_L = h_5(i) - h_6;
    
    % Efficiency
    eff(i) = (W_t_hp(i) + W_t_lp(i) - W_p) /Q_h(i);
end 
eff = eff *100;
x = linspace(150, 400, 6);
y = eff;

scatter(x, y, 'filled'), grid on
xlim ([100 450]); ylim([36 37])
xlabel('Reheat Presssure, (Psia)');
ylabel('Efficiency of Cycle');
title('Reheat Pressure vs Cycle Efficiency')

