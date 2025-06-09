close all, clc
%% #################### HW-5.3 ##################################
%% Givens
% the heat-flux profile form is
% q_ddt = q_ddt_0 *sin(pi(z +lam)/H_e);

q_ddt = 252500; % channel maximum heat flux, BTU/hr-ft^2
lam = 0.85;     % extrapolation distance, ft
H = 13.75;      %
p = 0.475/12;   % Pitch, ft
Do = 0.37/12;   % Outer clad dia, ft
Ro = Do/2;      % Outer radius, ft
D = 0.31/12;    % Pellet diameter, ft
R = D/2;        % Pellet radius, ft
t = 0.023/12;   % clad thickness, ft
k_c = 10;       % clad thermal conducitivity, BTU/hr-ft-F
H_G = 1000;     % Gap conductance, BTU/hr-ft^2-F
P = 2150;       % System Pressure, Psia
T_in = 575;     % Coolant inlet temp', F
T0 = 2090;      % Fuel centerline temp', F

H_e = H + 2*lam;
% = 15.45 ft

% Coolant properties at inlet temp' and system pressure can are
cp = 1.323;         % BTU/lb-F
k_cl = 0.321;       % BTU/hr-ft-F
mu = 0.2115;        % lb/ft

% Channel Area
A_x = p^2 - pi*Do^2/4;
% = 8.2017e-04 ft^2

% Prandtl number
Pr = cp * mu/k_cl;
% = 0.8717

% Wetted perimeter for flow
Pw = pi*Do;
% = 0.0969 ft

% Effective diameter
D_e = 4*A_x/Pw;
% = 0.339

%% Expressions
%%# Coolant Temperature 
% // write the equation herer

%%# Fuel Pellet Surface Temperature
% //

%%# Fluid Centerline Temperature
% //

%% Fuel Surface Temperature
% for fuel centerline temperature of T = 2090 F, interpolating from the coonducitivity integral
x11 = 2000;     x12 = 2100;
y11 = 5456.54;  y12 = 5611.14;
x1 = 2090;
y1 = intrp(x11, x12, y11, y12, x1);
% So,
int_k_T0 = y1;
% = 5.5957e+03 

int_k_TR = int_k_T0 - (q_ddt * Ro)/2;
% = 3.6493e+03

% Again looking into conductivity integral
x21 = 3561.45;  y21 = 1000;
x22 = 3791.87;  y22 = 1100;
x2 = int_k_TR;
y2 = intrp(x21, x22, y21, y22, x2);

% So, the fuel surface temperatur is-
T_R = y2;
% = 1.0381e+03 F

%% Coolant Temperature
syms G  % declaring varibale for mass flowrate

% So, the Coolant temperature  becomes, at H/2
T_cl = T_in + ( (q_ddt *Do *H)/(G *A_x *cp) );
T_cl = 9.8656e+07 *G + 575; %$ Redo  % Eq(1)

%% Convective HT Coefficient
% for Weisman Correlatin,
C = 0.042 *p/Do - 0.024;    % constant
% = 0.0299

% Reynolds Number
Re = G*D_e /mu;
Re = 0.1601 *G;  %$ Redo

% So, Convective HT Coefficent is 
h_cl = k_c/D_e *C *Re^0.8 *Pr^0.33;  % from Weisman Corr'
h_cl = 1.9498 *G^0.8; %$ Redo

% So, coolant temp' from other perspective
T_cl2 = T_R - q_ddt *( Ro/(R*H_G) + Ro/k_c *log(Ro/R) + 1/h_cl);
T_cl2 = 667.8923 - 1.2950e+05/G^0.8;    %$ Redo % Eq(2)

%% Mass Flux
% So, combining the values got in Eq(1) and Eq(2), we have-
eqn = T_cl == T_cl2;
% Here,
%   T_cl = 9.8656e+07 *G + 575
%   T_cl2 = 667.8923 - 1.2950e+05/G^0.8;

% Solving the above equation iteratively, we get mass flux
G2 = vpasolve(eqn, G, 1e6);
% = 1.315e06 lb/hr-ft^2

G3 = solfor(1.3e6, 1.5e6, 1e3);


%% ################################################################
% the interpolating function
function y = intrp(x1, x2, y1, y2, x)
    y = (x-x1)/(x1-x2) *(y1-y2) + y1;
end 

% function with for
function b = solfor(a,b, lim)
f = @(x)9.8656e+07 *x + 575 - 667.8923 + 1.2950e+05/x^0.8;
% a = accurate value
% b = initial guess
    for i=1:1:lim
       if f((a+b)/2) <0 % half-splitting
                a = (a+b)/2;
       else
                b = (a+b)/2;
       end
    end 
end 

% function with while
function b = solwhl(a, b, tol)
% a = accurate value
% b = initial guess
f = @(x)9.8656e+07 *x + 575 - 667.8923 + 1.2950e+05/x^0.8;
    while f(b) > tol % tolerance
        if f((a+b)/2)<0% % half-spliting
            a = (a+b)/2;
        else
            b = (a+b)/2;
        end
    end
end

