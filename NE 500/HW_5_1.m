clear all, close all, clc
%% ########## HW_5.1 ############################################
% The units are in SI unit
%% Givens
P = 1.5 *1e-2;  % pitch, m
D_f = 1.1 *1e-2; % fuel rod diameter, m
H = 3.6;        % fuel rod height, m
q_ddt = 90;     % heat flux, W/m^2
T_c_mx = 250;   % coolant max' temp', 0C
T_co_mx = 350;  % max' cladding surface temp, 0C' % design condition

    %% Material Properties
rh_w = 735.3;   % density, water, kg/m^3
rh_h = 0.865;   % density, helium
cp_w = 5.317 *1e3; % speficif heat, water, J/kg.K
cp_h = 5.23 *1e3; % speciffic heat, helium
mu_w = 0.955 *1e-4; % viscosity, water, kg/m.s
mu_h = 0.298 *1e-4; % viscosity, helium
k_w = 0.564;    % conductivity, water, W/m.K
k_h = 0.23;     % conductivity, helium


%% Calculation
A_f = P^2 - pi/4 *(D_f)^2; % flow area
% = 1.2997e-4 [m^2]
  
D_e = 4*A_f /(pi*D_f); % hydraulic diameter
% = 0.015 [m]

Pr_w = mu_w *cp_w /k_w; % Prandtl Number, water
% = 0.9

Pr_h = mu_h *cp_h /k_h; % Prantl number, heliuim
% = 0.678

syms m_w m_h    % creating variables

Re_w = m_w *D_e /(A_f * mu_w); % Reynolds no., water
Re_w = 1.2120e+06 *m_w;

Re_h = m_h *D_e /(A_f * mu_h); % Reynolds no., helium
Re_h = 3.8842e+06 *m_h;

% Now, nussselt number can be written as a function of mass flowrate and for this the appropriate
% one is Weissman Correlation

Nu_w = 0.023 *Re_w^0.8 * Pr_w^0.4 *(1.826 *(P/D_f) - 1.043 );   % Nusselt, Water
% and,
Nu_w = 2.3483e+03 *m_w^0.8;     % equation (I)

% Similarly
Nu_h = 0.023 *Re_h^0.8 * Pr_h^0.4 *(1.826 *(P/D_f) - 1.043 );   % Nusselt, Helium
Nu_h = 5.3215e+03 *m_h^0.8;     % equation (II)

% from the definition of Nusselt number
h_w = Nu_w *k_w /D_e; % convective HTC of water
h_w = 8.8041e+04 *m_w^0.8;      % equation (III)

h_h = Nu_h *k_h /D_e; % convective HTC of helium
h_h = 8.1360e+04 *m_h^0.8;      % equation (IV)

% As, h_w = convective coefficient
Del_T = q_ddt * (H*pi*D_f /(m_w*cp_w) + 1 /h_w);
    % Here, 
Del_T = T_co_mx - T_c_mx;   % The design temp' difference
    % = 100 0C

h_w = inv( Del_T/q_ddt - (H*pi*D_f) /(m_w*cp_w) );
h_w = inv(1.1111 - 2.3398e-05/m_w);
    a1 = 1.111; b1 = 2.3398e-05; c1 = b1/a1;
h_w = a1^-1 * (1 + c1 *m_w ); % binomial expansion
h_w = 0.9001 + 1.8956e-0 *m_w;   % equation (V)

% Again,
h_h = inv( Del_T/q_ddt -  (H*pi*D_f) /(m_h*cp_h) );
h_h = inv(1.111 - 2.3787e-05/m_h);
h_h = 0.9001 + 1.9271e-05 *m_h;   % equation (VI)

% Now, combining equation III, V, and IV, VI, we get-

% RHS of equation VII
eq_Rw = Nu_w *k_w /D_e;         % Eq R(VII) 
eq_Rw = 8.8041e+04 * m_w^0.8;   
% LHS of equation
eq_Lw = 0.9001 + 1.8956e-0 *m_w; % Eq L(VII)

% Again, RHS of equation VIII
eq_Rh = Nu_h *k_h /D_e;         % Eq R(VIII)
eq_Rh = 8.1360e+04 * m_h^0.8;
% LHS of equation VIII
eq_Lh = 0.9001 + 1.9271e-05 *m_h; % Eq L(VIII)

% combining these two for solving iteratively
eqn_w = eq_Rw == eq_Lw;
eqn_h = eq_Rh == eq_Lh;

% So, flowrate of water
m_dt_w = vpasolve(eqn_w, m_w, 0.3);
% = 0.289 kg/s

% flowrate of helium
m_dt_h = vpasolve(eqn_h, m_h, 0.3);
% = 0.299 kg/s

% Which shows that mass flowrate of helium is larger
return 
%% Doing Other Way
% for wwater
T1_w = q_ddt / Del_T;
T2_w = H*pi*D_f/cp_w;
T3_w = D_e /(2.3483e+03 *k_w); % Denom' of T3
            % constant comes from Nu_w
eqn_w1 = m_w == T1_w *(T2_w + T3_w *m_w^0.2);
m_dt_w = vpasolve (eqn_w1, m_w, 0.3);

% for helium
T1_h = q_ddt / Del_T;
T2_h = H*pi*D_f/cp_h;
T3_h = D_e /(5.3215e+03 *k_h);

eqn_h1 = m_h == T1_h *(T2_h + T3_h *m_h^0.2);
m_dt_h = vpasolve (eqn_h1, m_h, 0.3);










