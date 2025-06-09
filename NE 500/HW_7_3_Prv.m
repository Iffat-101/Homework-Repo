clear all, close all, clc
%% HW 7.3 ########## ########## ########## ########## 
Q = 4.2e9;
gam = 0.95;
a = 8;  lm_x = 2.30;    ae = a +2*lm_x;
b = 10; lm_y = 3.20;    be = b +2*lm_y;
c = 12; lm_z = 0.25;    ce = c +2*lm_z;
  
q_dt_mx = 14.5e3;   % KW/ft
q_2dt_mx = 475000;  % BTU/hr-ft^2
D = 0.3455/12;      % ft

%% a

F_nm = a*b*c;
F_dnm = 8*ae*be*ce/pi^3 *sin(pi*a /(2*ae)) *sin(pi*b /(2*be)) *sin(pi*c /(2*ce));
Fq = F_nm /F_dnm;

H = c; % height of core
n = Q*gam*Fq /(q_dt_mx *H);
% 48170

%% b
C = 4.412;      % conv. constant, W to BTU
D_o = Q*C *gam*Fq /(q_2dt_mx*pi*n*H);

%% c
q_3dt_0 = Q *gam*Fq /(n *pi*D^2/4 *H);
x = a/2;
q_3dt_0_rd = q_3dt_0 * cos(pi*x/(ae));

q = q_3dt_0_rd *pi*D^2/4 *2*ce/pi *sin(pi*c/(2*ce))

