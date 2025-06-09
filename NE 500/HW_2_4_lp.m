clear all, close all, clc
%% #################### HW-2.4(b) ##################################
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

%% Data for varying Tap Pressure
P_3 =  [35,      45,      55,      65,      75,      85,      95,      105,     115,     125];
T_3 =  [259.25,  274.41,  287.06,  297.96,  307.59,  316.24,  324.12,  331.36,  338.08,  344.35]; % no need
h_3f = [228.03,  243.5,   256.45,  267.65,  277.59,  286.55,  294.73,  302.27,  309.28,  315.85];
s_3f = [0.381,   0.402,   0.419,   0.434,   0.447,   0.459,   0.469,   0.479,   0.487,   0.496];
v_3f = [0.0171,  0.0172,  0.0173,  0.0174,  0.0175,  0.0176,  0.0177,  0.0178,  0.01785, 0.0179];
h_3g = [1167.19, 1172.16, 1176.11, 1179.36, 1182.12, 1184.49, 1186.55, 1188.37, 1189.99, 1191.45];
s_3g = [1.687,   1.667,   1.651,   1.638,   1.626,   1.616,   1.607,   1.599,   1.591,   1.585];

P_6 = P_3;
P_7 = P_3;
h_7 = h_3f; s_7 = s_3f; v_7 = v_3f;

%% Calculation
for i = 1:1:10
    % High pressure tap
        % steam quality at 3
        s_3(i) = s_2; 
        x_3(i) = (s_3(i) - s_3f(i)) /(s_3g(i) - s_3f(i));
    
        % enthalpy at 3
        h_3(i) = h_3f(i) + x_3(i) * (h_3g(i) - h_3f(i));
    
    % Turbine outlet    % non iterative
        % quality at 4
        s_4 = s_2; % isentropic expansion
        x_4 = (s_4 - s_5) /(s_4g - s_5);
    
        % enthalpy at 4
        h_4 = h_5 + x_4 *(h_4g - h_5);
    
    % Condenser Pump Work, Input
    W_cp_m3(i) = v_5 *(P_6(i) - P_5) /cf;
        % this is normalized by m3 = m1 - m2
        
        % enthalpy at 6
        h_6(i) = h_5 + W_cp_m3(i);
    
        % Feedwater Heater
        % by first law, we get mR = m1/m2
        mR(i) = (h_7(i) - h_6(i)) /(h_3(i) - h_6(i)); %# check it
        
        % Condenser pump work
        W_cp_m1(i) = (1 - mR(i)) * W_cp_m3(i);
            % this is normalized by m1
    
    % Feed Pump Work, Input
    W_fp_m1(i) = v_7(i) *(P_1 - P_7(i)) /cf;
    
        % enthalpy at 1
        h_1(i) = h_7(i) + W_fp_m1(i);
    
    % Boiler heat, input
    Q_h_m1(i) = h_2 - h_1(i);
    
    % Turbine Work, Output
    W_t_m1(i) = h_2 - h_3(i) * mR(i) - h_4 *(1-mR(i) );
    
    
    % Cycle Efficiency
        % Net work
        W_nt(i) = W_t_m1(i) - (W_cp_m1(i) + W_fp_m1(i) );
    
        % efficiency;
        eff(i) = W_nt(i) / Q_h_m1(i) *100;
end 

%% Plotting the graph
x = P_3; y = eff;
scatter(x, y), grid on
xlabel('Tap Pressure, (psia)');
ylabel('Efficiency of Cycle, (%)');
title('Tap Pressure vs Cycle Efficiency')

