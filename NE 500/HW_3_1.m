clear all, close all, clc
%% ########## HW-3.1 #######################################
%% Givens and Table Readings
  % Units: P = Psi, T = 0F 
  % h = BTU/lb, s = BTU/lb-F v = ft^3/lb

% Pressure Zones
P_a = 1240;     T_st_a = 571.43;
h_af = 577.39;  s_af = 0.776;   v_af = 0.0225;
h_ag = 1182.36; s_ag = 1.364;    

P_b = 940;      T_st_b = 537.19;
h_bf = 533.17;  s_bf = 0.734;
h_bg = 1194.78; s_bg = 1.398;
    
P_c = 580;      T_st_c = 482.61;
h_cf = 467.48;  s_cf = 0.667;   v_cf= 0.02;
h_cg = 1204.17; s_cg = 1.449;

P_d = 385;      T_st_d = 440.90;
h_df = 420.01;  s_df = 0.617;
h_dg = 1204.92; s_dg = 1.489;

P_e = 160;      T_st_e = 363.55;
h_ef = 336.09;  s_ef = 0.521;
h_eg = 1195.51; s_eg = 1.565;

P_h = 110;      T_st_h = 334.78;
h_hf = 305.84;  s_hf = 0.483;
h_hg = 1115.76; s_hg = 1.595;

P_i = 25;       T_st_i = 240.03;
h_if = 208.51;  s_if = 0.353;   v_if = 0.0169;
h_ig = 1160.55; s_ig = 1.714;
    
P_j = 2;        T_st_j = 126.02;
h_jf = 94.02;   s_jf = 0.175;   v_jf = 0.0162;
h_jg = 1115.76; s_jg = 1.919;   

% Pressure Zones, Primary Loop
P_p = 2350;     T_st_p = 659.08;

P_q = 2250;     T_st_q = 652.74;
                v_qf = 0.0269;
    
% others
eff_t = 0.95; % turbine efficiency
eff_p = 0.8;  % pump efficiency
Cf = 5.4;     % convrsn factor for IU

%% Calculations ####################
% point 1
P_1 = P_c; T_1 = 480;
h_1 = 464.44; s_1 = nan; v_1 = v_cf;

% point 2
P_2 = P_a;  x_2 = 0.95;
h_2 = 1152.11; s_2 = 1.334;

    %% High Pressure Turbine
% Point 3
P_3 = P_b;
h_3 = h_2;  % isenthalpic component
h_3f = h_bf; h_3g = h_bg;   h_3fg = h_3g - h_3f;
s_3f = s_bf; s_3g = s_bg;   s_3fg = s_3g - s_3f;

% Quality at 3
x_3 = (h_3 - h_3f) /h_3fg;

% Entropy at 3
s_3 = s_3f + x_3 *s_3fg;

    %% High Pressure Turbine Tap
% At point 4
s_4 = s_3; % constant entropy process
h_4f = h_df; h_4g = h_dg;   h_4fg = h_4g - h_4f;
s_4f = s_df; s_4g = s_dg;   s_4fg = s_4g - s_4f;

% Quality at 4
x_4 = (s_4 - s_4f) /s_4fg;

% Standard enthalpy at 4
h_4 = h_4f + x_4 *h_4fg;

% Actual enthalpy at 4
h_4a = h_3 - eff_t *(h_3 - h_4);

    %% High Pressure Turbine Outlet
% At point 5
s_5 = s_3;  % isentropic process
s_5f = s_ef; s_5g = s_eg; s_5fg = s_5g - s_5f;
h_5f = h_ef; h_5g = h_eg; h_5fg = h_5g - h_5f;

% Quality at 5
x_5 = (s_5 - s_5f) /s_5fg;

% Standard Enthalpy at 5
h_5 = h_5f + x_5 *h_5fg;

% Actual Enthalpy at 5
h_5a = h_3 - eff_t *(h_3 - h_5);

    %% Reheater
% At point 6
P_6 = P_e; T_6 = 560;
h_6 = 1304.63;  s_6 = 1.684;

% At point 7
P_7 = P_a; T_7 = T_st_a;
h_7 = h_af; s_7 = s_af;

    %% Low Pressure Turbine, Tap 1
% At point 8
P_8 = P_h; 
s_8 = s_6; % isentropic process
s_8f = s_hf; s_8g = s_hg; s_8fg = s_8g - s_8f;
h_8f = h_hf; h_8g = h_hg; h_8fg = h_8g - h_8f;

% Quality at 8
x_8 = (s_8 - s_8f) /s_8fg;

% Standard Enthalpy at 8
h_8 = h_8f + x_8 *h_8fg;
%     h_8 = 1265.44;   %$$ Check it later

% Actual Enthalpy at 8
h_8a = h_6 - eff_t *(h_6 - h_8);

    %% Low Pressure Turbine, Tap 2
% At point 9
P_9 = P_i;
s_9 = s_6;  % isentropic process
s_9f = s_if; s_9g = s_ig; s_9fg = s_9g - s_9f;
h_9f = h_if; h_9g = h_ig; h_9fg = h_9g - h_9f;

% Quality at 9
x_9 = (s_9 - s_9f) /s_9fg;

% Standard Enthalpy at 9
h_9 = h_9f + x_9 *h_9fg;

% Actual Enthalpy at 9
h_9a = h_6 - eff_t *(h_6 - h_9);

    %% Low Pressure Turbine Outlet
% At point 10
P_10 = P_j;
s_10 = s_6; % isentropic expansion
s_10f = s_jf; s_10g = s_jg; s_10fg = s_10g - s_10f;
h_10f = h_jf; h_10g = h_jg; h_10fg = h_10g - h_10f;

% Quality at 10
x_10 = (s_10 - s_10f) /s_10fg;

% Standard Enthalpy at 10
h_10 = h_10f + x_10 *h_10fg;

% Actual Enthalpy at 10;
h_10a = h_6 - eff_t *(h_6 - h_10);

    %% Condensate Pump
% At point 11
P_11 = P_j;
h_11 = h_jf; s_11 = s_jf; v_11 = v_jf;

% At point 12
P_12 = P_i;

% Condensate Pump Work, Standard
W_cp = v_11 *(P_12 - P_11) /Cf;

% Condensate Pump Work, Actual
W_cp_a = W_cp /eff_p;

% Enthalpy at 12
h_12 = h_11 + W_cp_a;

    %% Condensate Booster Pump
% At point 13
P_13 =  P_i;
h_13 = h_if;    v_13 =  v_if;

% At point 14
P_14 = P_c;

% Booster Pump Work, Standard
W_cpb = v_13 *(P_14 - P_13) /Cf;

% Booster Pump Work, Actual,
W_cpb_a = W_cpb /eff_p;

% Enthalpy at 14
h_14 = h_13 + W_cpb_a;

    %% Low Pressure Feed-Heater
% At point 15
P_15 = P_h; T_15 = T_st_h;
h_15 = h_hf;

% At point 16
P_16 = P_c; T_16 = 330;
h_16 = 301.65; s_16 = nan;

    %% High Pressure Feedheater
% At point 17
P_17 = P_d; T_17 = T_st_d;
h_17 = h_df;

    %% Feed Pumps
% At point 21
P_21 = P_a; v_21 = v_af;

% Feed Pump Work, Standard 
W_fp = v_1 *(P_21 - P_1) /Cf;

% Feed Pump Work, Actual
W_fp_a = W_fp /eff_p;

% Enthalpy at point 21
h_21 = h_1 + W_fp_a;

%% Primary Loop &&&&%%%%
% At point 18
P_18 = P_p; T_18 = 645;  % change
h_18 = 683.84;  % Table read

% At point 19
P_19 = P_q; T_19 = 575;
h_19 = 578.99;          v_19 =  0.0221;

% At point 20
P_20 = P_p; 

% RC Pump Work, Standard
W_rcp = v_qf *(P_20 - P_19) /Cf;  % Check

% RC Pump Work, Actual
W_rcp_a = W_rcp /eff_p;

% Enthalpy at point 20
h_20 = h_19 + W_rcp_a;

%% Mass Flow Rates ##########
% HP Heater
    % m4 /m1, mass ratio
mR4 = (h_1 - h_16) /(h_4 - h_17);

% Reheater
    % m3 /m1
mR3 =  (1 - mR4) *(h_6 - h_5a) /(h_2 - h_7 + h_6 - h_5a);

% Feed Pump
mR2 = 0.5;

% LP Heater
mR5 = ((h_14 - h_16) + mR3 *(h_7 - h_15) + mR4 *(h_17 - h_15)) /(h_15 - h_8a);
mR5 = abs(mR5);  % negative value ensues here

% Open Heater
mR6 = ( h_12 - h_13 + (mR3 + mR4 + mR5) *(h_15 - h_12)) /(h_12 - h_9);

% Steam Generators
mRp = mR2 *(h_2 - h_21) /(h_18 - h_19);

% across Core
mRc = 2 *mRp;

%% Cycle Efficiency
    % Boiler/ Reactor Core Heating
 Q_h_m1 = mRc *(h_18 -  h_20);  % normalized by m1 

     % HP Turbine Work, Output
 W_ta_hp_m1 = (1 - mR3) *h_3 - mR4 *h_4a - (1 - mR3 -mR4) *h_5a;

    % LP Turbine Work, Output
W_ta_lp_m1 = (1 - mR3 - mR4) *h_6 - mR5 *h_8a - mR6 *h_9a -(1 -mR3 -mR4 -mR5 -mR6) *h_10a;
    
    % Condensate Pump Work, Input
W_cp_a_m1 = (1 - mR3 - mR4 -mR5 -mR6) *W_cp_a;

    % Booster Pump Work, Input
W_cpb_a_m1 = W_cpb_a;

    % Feed Pump Work, Input
W_fp_a_m1 = W_fp_a;

    % Reactor Coolant Pump, Input
W_rcp_a_m1 = 2* mRp *W_rcp_a;

% Cycle Efficiency
W_net = ( W_ta_hp_m1 + W_ta_lp_m1) - (W_cp_a_m1 + W_cpb_a_m1 + W_fp_a_m1 + W_rcp_a_m1);
eff = W_net /Q_h_m1 *100;

