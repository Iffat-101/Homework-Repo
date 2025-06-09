clear all, close all, clc
%% #################### HW-1.5 ##################################
%% Givens and Table Readings
P_1 = 1150;
T_1 = [561.9,   596.9,   631.9,   666.9];
DT =  [0,       35,      70,      105]; % superheat
h_1 = [1186.44, 1227.67, 1260.62, 1289.39];
s_1 = [1.37,    1.41     1.443,   1.47];

P_2 = 2.5; T_2 = 134.4;

P_2g = 2.5; T_2g = 134.4;
h_2g = 1119.26; s_2g = 1.9;

P_3 = P_2g; T_3 = T_2g;
h_3 = 102.37; s_3 = 0.189; v_3 = 0.016; x_3 = 0;

P_4 = P_1;

P_4d = P_1; T_4d = T_1;

Cf = 5.4;

%% Calculation
for i = 1:1:4
    % Steam quality at 2
    s_2(i) = s_1(i); % isentropic process
    x_2(i) = (s_2(i) - s_3) /(s_2g - s_3);
    
    % Enthalpy at 2
    h_2(i) = h_3 + x_2(i) *(h_2g - h_3);
    
    % Turbine work, output
    W_t(i) = h_1(i) - h_2(i);
    
    % Condenser Cooling
    Q_l = h_2(i) - h_3;
    
    % Pump Work, input
        % Enthalpy at 4
    h_4 = h_3 + v_3 * (P_4 - P_3)/ Cf;
    W_p = h_4 - h_3;
    
    % Heating in Boiler
    Q_h = h_1(i) - h_4;
    
    % Net work
    W_nt(i) = W_t(i) - W_p;
    
    % efficiency
    eff(i) = W_nt(i) / Q_h;

    % Moisture content
    m_wt(i) = 1 - x_2(i);
end
eff = eff *100;

x = DT; y = eff;
scatter(x, y), grid on
xlabel('Superheat Temperature, (dT)');
ylabel('Efficiency of Cycle, (%)');
title('Boiler Superheat vs Cycle Efficiency')

