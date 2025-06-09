clear all, close all, clc
%% 8.3-b
rh_l = 999.21;      % density
rh_g = 1.205;       
g = 9.81;           % Grvity
C_D = 0.49;         % drag coefficient
D_xg = linspace(0.5,30,100); % bubble dia, in mm
D_g = D_xg *1e-3;      % ..in meter


v_r = sqrt( 4*(rh_l - rh_g) *D_g*g /(3*C_D *rh_l) );

figure(1)
plot(D_xg, v_r, 'k'), grid on
xlabel('Bubble Diameter (mm)');
ylabel('Rise Velocity (m/s)');
title('Bubble Dia vs Rise Velocity');

%% 8.3-c

mu_l = 1e-3;    % viscosity, liquid
gm = 7.27e-2;    % surface tension
Eo = (rh_l - rh_g)*g *D_g.^2 /gm;   % Eotvos  number
Re_b = rh_l .*D_g .*v_r /mu_l;     % Bubble reynolds


T1 = 16 ./Re_b .*(1 + 2 ./( 1 + 16 ./Re_b + 3.315./sqrt(Re_b) ));
T2 = 4 *Eo ./(Eo +9.5);
C_Dv = sqrt( T1.^2 + T2.^2);
v_rv = sqrt( 4*(rh_l - rh_g) *D_g*g ./(3*C_Dv *rh_l) );

figure(2)
plot(D_xg, v_rv, 'g'), grid on
xlabel('Bubble Diameter (mm)');
ylabel('Rise Velocity (m/s)');
title('Bubble Dia vs Rise Velocity at Variable CD')

%% 8.3-d
y1 = 1 ./Eo;
y2 = 1 ./Re_b;

figure(3)
plot(D_xg, y1, 'b'), grid on
xlabel('Bubble Diameter (mm)');
ylabel('1/Eo', 'Interpreter','latex');
title('1/Eo vs Bubble Diameter')

figure(4)
plot(D_xg, y2, 'b'), grid on
xlabel('Bubble Diameter (mm)');
ylabel('$1/Re_b $', 'Interpreter','latex');
title('1/Re vs Bubble Diameter')

