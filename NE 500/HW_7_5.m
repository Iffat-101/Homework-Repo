clear all, close all, clc
%% HW 7.5 ########## ########## ########## ########## 

m_Pu = 238.0495;    % mass, Plutonium, amu (gm/mole)
m_U = 234.0409;     % mass, Uranium, amu 
m_He = 4.0026;      % Mass, Helium
m_O = 15.994;       % Oxygen

sg = 5.67e-8;       % Boltzman constant
k = 3.1;            % Conductivity of PuO2
T_mlt = 3000;
eff = 0.08;         % electrical's efficiency

%% Calc
Dm = m_Pu - (m_U + m_He);    % mass defect
Dm = Dm *931;       % Mev

rh_PuO2 = 11.5;  % Density, g/cm^3
Na = 6.023e23;  % Avogadro's constant
m_PuO2 = m_Pu + 2*m_O; 

N_PuO2 = rh_PuO2 *Na /m_PuO2; % Number density

t_h = 86 *365*24*3600;  % half life, sec
lam = log(2) /t_h;

Rt = lam *N_PuO2;        % decay rate, #/cm^3.s

q_3dt = Dm *Rt;          % decay-energy stored, MeV/cm^3.s
C_MeV_J = 1.602e-13;    % MeV to Joules
C_cm3_m3 = 1e-6;        % cc to m^3
q_3dt2 = q_3dt * C_MeV_J/C_cm3_m3;    % Watt/m^3

T_lim = 0.85*T_mlt;

syms x
eqn = (q_3dt2 *x /(3*sg) )^0.25 + q_3dt2 *x^2 /(6*k) == T_lim;
R = vpasolve(eqn, x, 0.007);    % 0.007 initial parameter
R = double(R);

V = 4/3*pi *R^3;       % volume of sphere
Q_el = q_3dt2 *V *eff;
