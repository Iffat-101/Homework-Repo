clear all, close all, clc
%% 8.4-a
rh_l = 999.21;      % density
rh_g = 1.205;       
g = 9.81;           % Grvity
C_D = 0.49;         % drag coefficient
D_xg = [1, 1.5, 2];       % bubble dia, in mm
D_g = D_xg *1e-3;      % ..in meter
alf = [0.03, 0.01, 0.01];

mu_cl = 1.05e-6;    % viscosity, continuous liquid

v_r = sqrt( 4*(rh_l - rh_g) *D_g *g /(3*C_D *rh_l) );

mu_2p = 0.6 .*D_g .*alf .*v_r;
mu_2p_n = mu_2p /mu_cl; % normalized


%% 8.4-b
alf_e = 0.05;        % Equivalent VF
D_ge = 3/5*D_g(1) + 1/5*D_g(2) + 1/5*D_g(3);

v_re = sqrt( 4*(rh_l - rh_g) *D_ge *g /(3*C_D *rh_l) );
mu_2pe = 0.6 .*D_ge .*alf_e .*v_re;
mu_2pe_n = mu_2pe /mu_cl;
