clear all, close all, clc
%% #################### HW-2.2 ##################################
% Units: P = KPa, T = K 
% h = J/kg, s = J/kg-k, v = m^3/kg

% at point 1
P_1 = 6880 ; 
h_1 = 2774.11; 

% at point 2
P_2 = 3;%% (kPa) 0.03 (bar)
s_1 = 5.82;
s_f = 0.35;
s_g = 8.57;
s_fg = s_g - s_f;

% direct from the tablex2s
x_2 = (s_1 - s_f)/s_fg;
h_f = 100.99;
h_g = 2544.87;
h_fg = h_g-h_f;
h_2 = h_f + x_2* h_fg;

%h_1=2785.9;
W_ts = (h_1 - h_2); % HP turbine work
W_ta = 0.9 *W_ts;
h_2a = h_1 - W_ta;
x_2a = (h_2a - h_f) /h_fg;

% At point 3 
P_3 = 3;      
h_3 = 100.99; 

h_4s = 192.6;
% W_ps = (h_3 - h_4);

eff_P=0.8;
v_3 = 0.001;

P_4 = 6880;
W_ps = v_3 *(P_4- P_3);
W_pa = W_ps/eff_P;
h_4a = h_3 + W_pa; 

 % Boiler heating
Q_h = h_1 - h_4a; %%

% Efficiencies
eff_wo_m = (W_ta - W_pa)/Q_h;
eff_w_m = 0.3611;
% Efficiency increment
eff_inc = (eff_w_m - eff_wo_m);
 