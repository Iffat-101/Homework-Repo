clear all, close all, clc
%% ########## HW-6.3 ############################################
% Here, all the units are in English Units
%% Given Values
% Tube Side
Do = 0.90/12;       % Tube outer diameter, ft
t = 0.045/12;       % Tube thickness, ft
Di = Do -2*t;        % Tube inner diameter, ft
 % = 0.0675 ft
k_w = 209;          % Tube thermal conductivity, BTU/hr-ft.F
n_t = 440;          % no. of tubes
m_t = 3.7e+05;      % Mass flowrate, lb/hr
T_h1 = 120;         % Inlet temp', F
T_h2 = nan;         % Exit Temp', F,    TBD
% Assuming sat' temp' at T_h1
cp_t = 0.998;       % Specific heat, BTU/lb.f
mu_t = 1.348;       % Viscosity, lb/ft.h
k_t = 0.369;        % Thermal conductivity, fluid's, BTU/hr.ft.F

% Shell Side
s = 1.25/12;        % Tube pitch, ft
m_s = 1.90e+05;     % Mass flowarate, lb/hr
T_c2 = 70;          % Inlet temp', F
T_c1 = nan;         % Exit temp, F      TBD
% Assuming sat' temp' at T_c2
cp_s = 0.999;       % Specific heat, BTU/lb.F
mu_s = 2.358;       % Viscosity, lb/ft.h
k_s = 0.347;        % Thermal conductivity, BTU/hr.ft.F

% Lattice type, square

Q_dt = 1.50e+06;    % Capacity of exchanger, MW
Q_dt = Q_dt *3.412; % in BTU/hr
% = 5118000 BTU/hr

%% Finding Missing Values
% As we know
% Q_dt = m_t *cp_t (T_h1 - T_h2) 
T_h2 = T_h1 - Q_dt /(m_t *cp_t);
% = 106.1398 F

% Similarly,
T_c1 = T_c2 + Q_dt /(m_s *cp_s);
% = 96.9638 F

% Tube side avg. temp
T_h_av = (T_h1 + T_h2) /2;
% = 113.0699 F

% Shell side avg. temp'
T_c_av = (T_c1 + T_c2) /2;
% = 83.4819 F

%% LMTD
% //eqn
  Dl_Tm_Nm = ( T_h2 - T_c2) - (T_h1 - T_c1);        % Numerator
  Dl_Tm_Dn = log ( (T_h2 - T_c2) /(T_h1 - T_c1) );  % Denominator
Dl_Tm = Dl_Tm_Nm/ Dl_Tm_Dn;
% = 29.0979 F

%% Finding UA
UA = Q_dt /Dl_Tm;
% = 1.7589e+05  BTU/hr.F

%% Convective HTC, Tube Side
% Flow area
A_x_t = n_t *pi *Di^2 /4;
% = 1.5745 ft^2

% Mass flux
G_t = m_t/ A_x_t;
% = 2.3499e+05  lb/ft^2.hr

% Reynolds number
Re_t = G_t *Di /mu_t;
% = 1.1767e+04

% Prandtl number
Pr_t = cp_t *mu_t /k_t;
% = 3.6458

% Nusselt Number
Nu_t = 0.023 *Re_t^0.8 *Pr_t^0.3; % for cooling fluid
% = 61.2066

% Convective HTC
h_t = k_t /Di *Nu_t;
% = 334.5962 BTU/hr.ft^2.F

%% Convective HTC, Shell Side
% As, it is a square lattice
% Flow Area
A_x_s = s^2 - pi *Do^2 /4;
% = 0.0064 ft^2

% Wetted Perimeter
P_w = pi *Do;
% = 0.2356 ft

% Effective Diameter
D_e = 4 *A_x_s/P_w;
% = 0.1092 ft^2

% Mass flux
G_s = m_s /(n_t *A_x_s);
% = 6.7127e+04 lb/hr.ft^2

% Reynlods number
Re_s = G_s *D_e /mu_s;
% = 3.1089e+03

% Prandtl number
Pr_s = cp_s *mu_s /k_s;
% = 6.7886

% Nusselt Number
    % Weissman Correlation, Rectangular Channel
Nu_s = ( 0.042 *s/Do - 0.024 ) *Re_s^0.8 *Pr_s^0.33;
% = 40.2057

% Convective HTC
h_s = k_s /D_e *Nu_s;
% = 127.7517 BTU/hr.ft^2.F

%% Heat Exchanger Length
L = UA/n_t *( 1/(h_s *pi*Do) + 1 /(2 *pi *k_w) *log(Do/Di) + 1/(h_t *pi*Di) );
% = 18.9463 ft



