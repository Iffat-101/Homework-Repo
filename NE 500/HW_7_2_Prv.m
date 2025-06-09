clear all, close all, clc
%% HW 7 ########## ########## ########## ########## 
H = 150/12;         % core height 
lm_z = 3.25/12;     % extrapolation
lm_r = 24.2/12;
alf = 3.965;        % ft^-1
D = 0.33/12;        % pellet diameter, ft
D_r = 0.385/12;     % rod diameter, ft
D_c = 144/12;       % effective core dia, ft
Q  = 3.95e9;        % Core thermal output, W
q_3dt_c = 21.5e6;   % core volumetric heat gen', W/ft^3
q_3dt_c2 = q_3dt_c *3.412; % BTU/h
gam = 0.95;         % fuel heat fraction
n = 56900;          % no. of fuel rod
Fq = 2.35;          % total power peaking factor
Fz = 1.45;          % axial power peaking factor

%% a
R = D/2;
R_r = D_r/2;        % fuel rod diameter
R_c = D_c/2;        % core radius

He = H + 2*lm_z;    % extrapolated core/rod height
Re_c = R_c + lm_r; % extrapolated core radius;

J1 = besselj(1, alf*R);

C1 = 4*R*He /alf *J1 *sin( pi*H /(2*He));
q = q_3dt_c * C1;

%% b
q_dt_mx = q_3dt_c *2*pi *R/alf *J1;

q_dt_mx2 = Q *gam*Fq /(n*H); % not needed though

%% c
q_2dt_mx = q_dt_mx /(pi*D_r);

%% d
q_2dt_av = q_2dt_mx /Fq;
q_dt_av = q_2dt_av *pi *D_r;

q_3dt_av = q_dt_av *H /C1;

syms x
q_3dt_R = q_3dt_av /q_3dt_c;    % ratio of volumetric heats
Cb = 2.405/Re_c;                  % const. for Bessel

eqn = q_3dt_R == besselj(0, Cb*x);
r_av = vpasolve(eqn, x, 5.0);   % guess value, 5.0