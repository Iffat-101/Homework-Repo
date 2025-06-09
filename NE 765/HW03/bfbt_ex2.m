close all, clear all, clc
%% Input Values
P_o = 7.15;         % outlet pressure; MPa
T_av = 284.9;       % Avg temp; dC
m_fl = 20.3;        % mass inflow; t/hr
wd = 132.5;         % assembly width; mm
rd = 8.0;           % assembly radius; mm
cf = 1.0;           % friction coeff
ksp = 1.2;          % spacer-grid p-loss factor

%% Varying P
% P_o1 = 7.20;        % outlet pressure; MPa
% P_o2 = 7.10;
% n = 10;
% P_o = linspace(P_o2, P_o1, n); 
% 
%     % Model DP
% n = length(P_o);
% for i = 1:1:n
%     DP(i) = bfbt_T(P_o(i), T_av, m_fl, wd, rd, cf, ksp);
% end
% DP

%% Varying T
T_av1 = 287.0;       % Avg temp; dC
T_av2 = 250.0;       % 
n = 10;
T_av = linspace(T_av2, T_av1, n); 

    % Model DP
n = length(T_av);
for i = 1:1:n
    DP(i) = bfbt_T(P_o, T_av(i), m_fl, wd, rd, cf, ksp);
end
DP

%% Variying FLowrate
% m_fl1 = 21.0;       % mass inflow; t/hr
% m_fl2 = 20.0;       %
% n = 10;
% m_fl = linspace(m_fl1, m_fl2, n); 
%     
%     % Model DP
% n = length(m_fl);
% for i = 1:1:n
%     DP(i) = bfbt_T(P_o, T_av, m_fl(i), wd, rd, cf, ksp);
% end
% DP;

%% Model DP
% DP


