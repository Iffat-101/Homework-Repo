clear all, close all, clc
%% ########## HW_5.2 ############################################
%% Givens
d_rn = 2.5 *1e-2;   % round tube dia, m
d_sq = 2.5 *1e-2;   % square tube side, m
H = 3.5;            % tube length, m

%% Matter Props
rh = 914;           % density, kg/m^3
mu = 8.966 *1e-5;   % viscosity, Pa.s
cp = 5.49 *1e3;    % specific heat, J/kg-C
k = 0.565;          % conductivity, W/m-C
m_dt = 2.95 *1e-3;  % mass flowrate, kg/s

%% Calculation
%% # 01
    %# for Round Tube
A_rn = pi*d_rn^2 /4;    % flow area, round pipe, m^2
% = 4.9087e-04 m^2
P_w_rn = pi *d_rn;     % Wetted perimeter, m
% = 0.0785 m

D_e_rn = 4 *A_rn/P_w_rn;   % Effective diameter, m
% = 0.025 m

% So, Reynolds number
Re_rn = m_dt*D_e_rn /(A_rn*mu);
% = 1.6757e+03;
    % which is less than 2300, so, laminar flow

    %# for Square Tube
A_sq = d_sq^2;      % flow area, square pipe, m^2
% =  6.2500e-04 m^2

P_w_sq = 4 *d_sq;
% = 0.1 m

D_e_sq = 4 *A_sq/P_w_sq;
% = 0.0250 m

Re_sq = m_dt*D_e_sq /(A_sq *mu);
% =  1.3161e+03
    % which is less than 3, so laminar flow
 
%% # 02
% // insert a pic here

% Fully developed flow legnth
L = 0.05 *Re_rn * d_rn;
% = 2.0946 m
% so, the flow is not fully developed

%% # 03
% for fully developed region
Pr = cp *mu/k;
% = 0.8712

% As it is assumed as constant temp' surface condition
Nu_rn = 3.66;   % for round tube
Nu_sq = 2.98;   % for square tube

% As we know-
% Nu = h *D_e /k; 

h_rn = Nu_rn *k /D_e_rn;    % Heat transfer coeff. for round tube, W/m^2.k
% = 82.7160 W/m^2.k

h_sq = Nu_sq *k /D_e_sq;    % HTC for square tube
% = 67.3480 W/m^2.k

% Here, as teh round tube has better HTC so, it is good

%% # 04
m_dt3 = m_dt *3;    % tripled flowrate
% = 0.0089 kg/s

% So, the reynolds no. are
Re_rn3 = m_dt3 *D_e_rn /(A_rn *mu);
% = 5.0271e+03

Re_sq3 = m_dt3 *D_e_sq /(A_sq *mu);
% = 3.9482e+03

% So, in both cases, trippling the flowrate makes them into turbulent regime
% So, the n usselt number would be as follows
Nu_rn3 = 0.023 *Re_rn3^0.8 * Pr^0.4;
% = 19.8989
Nu_sq3 = 0.023 *Re_sq3^0.8 * Pr^0.4;
% = 16.4021

% then HTCs are
h_rn3 = Nu_rn3 *k /D_e_rn;
% = 449.7147 W/m^2.k

h_sq3 = Nu_sq3 *k /D_e_sq;
% = 370.6884 W/m^2.k

% Here, the round tube is also better than the square tube