clear all, close all, clc
%% ########## HW-6.2 ############################################
% Here, all the units are in English Units
%% Given Values
    %% Fluid Prop's
% //pic
% Primary Side, Within tube
m_1 = 79e+06;        % Mass flowrate, lb/hr
P_h = 2150;         % Pressure, Psia
T_h1 = 630;         % Hot leg temp', F
h_h1 = 659.77;      % Enthalpy at inlet, BTU/lb
h_h2 = nan;         % Enthalpy at exit, BTU/lb, TBD
T_h2 = nan;         % Primary exit temp', F,    TBD
cp_h = 1.69;        % Primary specific heat, BTU/lb.F

% Secondary Side, within shell
m_2 = 5.9e+06;       % Mass flow-rate, lb/hr
P_c = 930;          % Pressure, Psi
T_c1 = 575;         % Exit temp', F
T_c2 = 535.92;      % Inlet temp', Saturated, F
h_c2 = 1195.13;     % Enthalply at 'inlet', BT/lb
h_c1 = 1235.54;     % Enthalpy at exit, BTU/lb
 
T_c_av = (T_c2 + T_c1)/2;
% = 555.4600 F      % Average temp. of flow, F
% Now, evaluating values at T_c_av
rho_c = 1.96;       % Vapor density, lb/ft^3
cp_c = 1.02;        % Vapor specific heat; BTU/lb.F
mu_c = 0.046;       % Vapor viscosity, lb/ft.hr      
k_c = 0.034;        % Vapor thermal conductivity, BTU/hr-ft.F 

    %% Geometric
k_w = 10.7;         % Thermal conductivity of tube, BU/hr-ft.F
n_t = 15670;        % no. of tubes
s = 1.3 /12;        % pitch, ft
Do = 0.635 /12;     % Outer diameter of tube, ft
Di = 0.545 /12;     % Inner diameter of tube, ft
% Lattice is triangular

%% Evaluating Hot's Exit
% We get that
% Q = m2 *(h_c1 - h_c2) = m1 *(h_h1 - h_h2)

% So,
h_h2 = h_h1 - m_2 *(h_c1 - h_c2) /m_1;
% = 656.7520 BTU/lb

% As we know, Dl_T *cp_h = Dl_h
% Hot leg exit temp'
T_h2 = T_h1 - (h_h1 - h_h2) /cp_h;
% = 628.2142    F

% Using this value, we can further evaluate some flow props within primary side
T_h_av = (T_h1 + T_h2)/2;    % Primary average temp'
% = 629.1071    F

mu_h = 0.181;       % Primary viscosity, lb/ft.hr
k_h = 0.282;        % Primary thermal conductivity, BTU/h.ft.F

%% Finding LMTD
% //eqn
  Dl_Tm_Nm = ( T_h2 - T_c2) - (T_h1 - T_c1);        % Numerator
  Dl_Tm_Dn = log ( (T_h2 - T_c2) /(T_h1 - T_c1) );  % Denominator
Dl_Tm = Dl_Tm_Nm/ Dl_Tm_Dn;
% = 72.0455 F

%% Finding Q
% Heat Exchanged
Q_dt = m_1 *(h_h1 - h_h2);
% = 2.3842e+08  BTU/hr

%% Finding UA
% the product of Overall HTC (U) and area A is
UA = Q_dt /Dl_Tm;
% = 3.3093e+06  BTU/hr.F

%% Convective HTC, Tube Side 
% Flow area
A_x_t = n_t *pi*Di^2 /4;
% = 25.3857 ft^2

% Mass flux
G_t = m_1/ A_x_t;
% = 3.1120e+06  lb/ft^2.hr

% Reynolds Number
Re_t = G_t *Di /mu_h;
% = 7.8086e+05

% Prandtl Number
Pr_t = cp_h *mu_h /k_h;
% = 1.0847

% Nusselt Number
Nu_t = 0.023 *Re_t^0.8 *Pr_t^0.3; % for cooling fluid
% = 1.2201e+03

% HTC Tube Side   % don't mess-up between enthalpy and HTC
h_t = k_h /Di *Nu_t;
% = 7.5756e+03  BTU/hr.ft^2.F 

%% Convective HTC, Shell Side
% As this is a triangular lattice, So,
% Flow Area
A_x_s = 3^0.5/4 *s^2 - pi *Do^2 /8;
% = 0.0040 ft^2

% Wetted Perimeter
P_w_s = pi *Do/2;
% = 0.0831 ft

% Effective Diameter
D_e = 4 *A_x_s /P_w_s;
% = 0.1916 ft

% Mass flux
    % As for triangular lattice there're 2 flow channels per tube
G_s = m_2 /(2 *A_x_s) /n_t;
% = 4.7274e+04  lb/ft^2.hr

% Reynolds number
Re_s = G_s *D_e /mu_c;
% = 1.9694e+05

% Prandtl Number
Pr_s = cp_c *mu_c /k_c;
% = 1.3800

% Nusselt  Number
    % From Weissman Correlation, Triangular Setup
Nu_s = ( 0.026 *s/Do - 0.006 ) *Re_s^0.8 *Pr_s^0.33;
% = 903.3106

% HTC, Shell side
h_s = k_c /D_e *Nu_s;
% = 160.2650    BTU/hr.ft^2.F

%% Superheating Length
% Length of Superheating Area
L = UA/n_t *( 1/(h_s *pi*Do) + 1 /(2 *pi *k_w) *log(Do/Di) + 1/(h_t *pi*Di) );
% = 8.6020 ft

