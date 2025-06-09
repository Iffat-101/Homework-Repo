close all, clear all, clc
%% Input Values
P_o = 7.15;         % outlet pressure; MPa
T_av = 284.9;       % Avg temp; dC
m_flow = 20.30;     % mass inflow; t/hr

%% Inputs, minor
c_rh = 1;       % density adjustment factor
cv = 1;         % viscosity adjustment f
asm_wd = 132.5;     % assembly width; mm
asm_rd = 8.0;       % assembly radius; mm

htr = 6.15;         % g-parameters, height; mm
wtr = 17.0;         % g-parameter, width, mm
cf = 1.0;         % friction coeff
ksp = 1.2;          % spacer-grid p-loss factor
dZ = 3708e-03;      % domain length; m

%% Water props
T = T_av;           % temperature; dC
P = P_o*10;         % pressure, bar

rho = XSteam('rho_pT',P,T);     % density
rho = c_rh *rho;

mu = XSteam('my_pT',P,T);
mu = cv *mu;

%% Other Parameters
kgS_tHr = 1e3/60^2;
m_flow = m_flow *kgS_tHr;   % mass flow; kg/s

n = 100;     % no. of division
Cz = linspace(0, dZ, n+1);
dl = Cz(2) - Cz(1); % division
dz = dl*ones(1,n);

  asm_d =  2*asm_rd;      % assembly dia; mm
  asm_df = asm_wd - asm_d;  
A_BWR = asm_df^2 + 2*asm_d*asm_df + pi*asm_rd^2;    % mm^2
A = A_BWR - 60*pi*htr^2 - pi*wtr^2;       % effective area; mm^2

rg_m = 1e-3;        % conversion factor
A = A*(rg_m)^2;     % Area; m^2

Pw = 60*pi*12.3 + pi*34 + 4*asm_df + 2*pi*asm_rd;    % wetted p; mm
Pw = Pw*rg_m;   % wetted perimeter; m

Gm = m_flow /A;      % flux; kg/s.m^2
De = 4*A/Pw;        % effective diameter; m

v = m_flow /(rho *A);   % velocity; m/s
Re = Gm *De /mu;
f = cf *0.184 *Re^-0.2; % Reynold's friction-loss; -

dP_frm = ksp *rho*v^2/2;    % P-drop; Pa
       % P-loss due to spacer grid 

K_rg = 1e-3;                % conv f, rg-K

%% Looper
C1 = f*(Gm^2) /(2*De*rho);    % evaluating a const.

for i=1:1:n
    dP(i) = C1*dz(i);
    dP(i) = dP(i) *K_rg;        % KPa
end

dP;
CdP = cumsum(dP);


