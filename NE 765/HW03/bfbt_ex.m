close all, clear all, clc
%% Input Values, 1
P_o1 = 7.15;        % outlet pressure; MPa
T_av1 = 284.9;      % Avg temp; dC
m_fl1 = 20.3;       % mass inflow; t/hr
wd = 132.5;         % assembly width; mm
rd = 8.0;           % assembly radius; mm
cf = 1.0;           % friction coeff
ksp = 1.2;          % spacer-grid p-loss factor

y1 = bfbt_T(P_o1, T_av1, m_fl1, wd, rd, cf, ksp);
Y1 = bfbt_ax(P_o1, T_av1, m_fl1, wd, rd, cf, ksp);

fprintf('Total axial pressure drop [KPa]; case 1: %d \n', y1);
fprintf('Axial pressure-drop profile; case 1:' );Y1

%% Input Values, 2
P_o2 = 7.16;        % outlet pressure; MPa
T_av2 = 285.1;      % Avg temp; dC
m_fl2 = 24.9;       % mass inflow; t/hr
wd = 132.5;         % assembly width; mm
rd = 8.0;           % assembly radius; mm
cf = 1.0;           % friction coeff
ksp = 1.2;          % spacer-grid p-loss factor

y2 = bfbt_T(P_o2, T_av2, m_fl2, wd, rd, cf, ksp);
Y2 = bfbt_ax(P_o2, T_av2, m_fl2, wd, rd, cf, ksp);

fprintf('Total axial pressure drop [KPa]; case 2: %d \n', y2);
fprintf('Axial pressure-drop profile; case 2:' );Y2
