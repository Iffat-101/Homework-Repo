clear all,
close all, clc

%% P = 13 bar
rho_l = 874.27;
rho_g = 6.61 ;
Dl_rho = rho_l - rho_g;
sg = 39.58 *10^-3 ;
x = (Dl_rho)/(sg * 25906.6);

%% P = 20 bar
rho_l2 = 849.79  ;
rho_g2 = 10.04  ;
Dl_rho2 = rho_l2 - rho_g2;

sg2 = 34.83 *10^-3 ;
x2 = (Dl_rho2) /(sg2 * 25906.6);

% %% P = 25 bar
% rho_l2 = 835.12  ;
% rho_g2 = 112.58 ;
% Dl_rho2 = rho_l2 - rho_g2;
% 
% sg2 = 32.15 *10^-3 ;
% x2 = (Dl_rho2) /(sg2 * 25906.6)