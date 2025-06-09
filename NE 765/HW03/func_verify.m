close all, clear all, clc
%% Plot Function
x = linspace(0,2, 101);
y1 = x.^3;
y2 = 2.*x;

% plot(x,y1, x, y2), grid on

%% Analytical Definite Integral
syms x
y1 = x.^3;
y2 = 2.*x;
Y = y2 - y1;

eqn = y2 - y1 == 0;
S = solve(eqn,x);
S = S(2)

A = int(Y,0, S)

%% Numerical Part
x = linspace(0,2, 1001);
  L = length(x);
dx = x(2) - x(1);

y1 = x.^3;
y2 = 2.*x;

%%  The functions comes here
Df1 = y2 - y1;
  sum1 = 0;
for i=1:1:L-1  % trapezoidal method reduces 1 iteration
    if Df1(i) > 0        % fetching only positive values
        Y(i) = (Df1(i) + Df1(i+1))/2 *dx;  % trapezoidal method
        sum1 = Y(i) + sum1;
    else
        continue
    end
end
sum1 = round(sum1, 3)      % rounding up

%% #### 
function y = func(x)
    y1 = x.^3;
    y2 = 2.*x;
    y = y2 - y1;
end
