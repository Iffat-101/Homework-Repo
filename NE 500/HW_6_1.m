clear all, close all, clc
%% ########## HW-6.1 ############################################
% Here, all the units are in English Units
%% Design Variables,    Input
Do = 1.00 /12;      % Tube outer diameter, ft
n_t = 1000;         % number of tubes

%% Given Values
% //pic
% The shell side temperature is assumed as saturation temp'
% So,
P_h1 = 170;         % Shell-side pressure, Psia
T_h1 = 368.43;      % Shell-side temp', F   ; Table-read
T_h2 = T_h1;        % Shell-side outlet temp'

P_c = 285;          % Tube-side pressure, psia
T_c1 = 305;         % Tube-side inlet temp', F
T_c2 = 367;         % Tube-side outlet temp', F
T_c_av = (T_c1 + T_c2)/2; % Tube-side avg. temp'
% = 336 F
cp_t = 1.042;       % Specific heat at tube-side, BTU/lb.F
rho_t = 56.15;      % Fluid density at tube-side, lb/ft^3
mu_t = 0.389;       % Dynamic viscosity at tube side, lb/ft.h
k_t = 0.393;        % Thermal conductivity, tube fluid, BTU/hr-ft.F

% Do = 0.75 /12;      % Tube outer diameter, ft
t = 0.045/12;       % Tube thickness, ft
Di = Do - 2*t;      % Tube inner diameter, ft
% = 0.0758 ft

% n_t = 1000;         % number of tubes
m_t = 11.8e6;       % Tube-side mass flowrate, lbm/hr
k_w = 10;           % Tube-wall thermal conductivity, BTU/hr.ft.F

%% Finding LMTD
% //eqn
  Dl_Tm_Nm = ( T_h2 - T_c2) - (T_h1 - T_c1);        % Numerator
  Dl_Tm_Dn = log ( (T_h2 - T_c2) /(T_h1 - T_c1) ); % Denominator
Dl_Tm = Dl_Tm_Nm/ Dl_Tm_Dn;
% = 16.3491 F

%% Finding Q
% Heat transfer in tube-side
Q_dt = m_t *cp_t *(T_c2 - T_c1);
% = 762327200   BTU/hr

%% Finding UA
% the product of Overall HTC (U) and area A is
UA = Q_dt /Dl_Tm;
% = 4.6628e+07  BTU/hr.F

%% Convective HTC, Tube Side
% Flow area
A_x = n_t *pi *Di^2 /4;
% = 4.5166 ft^2

% Mass flux
G_t = m_t/ A_x;
% = 2.6126e+06  lb/ft^2.hr

% Flow velocity
v_t = G_t /rho_t;
% = 4.6529e+04  ft/hr

% Reynolds Number
Re_t = G_t*Di /mu_t;
% = 5.0931e+05

% Prantl Number
Pr_t = cp_t *mu_t /k_t;  % check
% = 1.0314 

% Nusselt Number
Nu_t = 0.023 *Re_t^0.8 *Pr_t^0.4;
% = 856.4138

% So, the HTC
hi = Nu_t *k_t /Di;
% = 4.4383e+03  BTU/ft^2.F.hr

%% Heat Transfer Length
% Now, we know
% //eqn
% Which means
% //eqn
ho = Inf;      % Initializing some value for shell-side convective HTC
L = UA/n_t *( 1/(ho *pi *Do) + 1 /(2 *pi *k_w) *log(Do/Di) + 1/(hi *pi *Di) );

% Considering negligible thermal resistance across tube-wall to shell side
L = UA/n_t *( 1 /(2 *pi *k_w) *log(Do/Di) + 1/(hi *pi *Di) )
% = 114.0872 ft

%% Heat Transfer Area
% Average diameter of tube
D_av = (Do + Di) /2;
% = 0.0796 ft

A = pi *D_av *L *n_t
% = 2.8524e+04  ft^2

%% Pressure Drop
% For the Reynolds number found, the flow is turbulent, so, 
% Friction Factor
f_D = 0.0055 *( 1 + 100 /Re_t^0.333);
% = 0.0124

% Pressure Drop can be found using Darcy-Weisbach Equation
% Pressure Drop
Dl_P = L * f_D *rho_t/2 *v_t^2/Di;
% = 1.1355e+12  lb/ft.hr^2

% Converting to psia
Dl_P2 = Dl_P * (1/32.2) * (1/60^2)^2 *(1/12^2)
             % lb->lbf   hr->s       ft^2->in^2
% = 18.8949 Psia


%% Creating a Table
% //table

