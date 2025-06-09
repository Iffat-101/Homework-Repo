clear all, close all, clc

%% Experimental Conditons
DP_e = [3.12, 3.00, 2.88, 3.18, 3.17, 2.92, 3.23, 2.84, 3.19, 3.15];

b0 = max(DP_e);
a0 = min(DP_e);

mu = (b0 + a0)/2;   % taking the interval
sg = 1.0;           % assumed STD, %## input

%% Grid
N = 100;            % no of division;
b = b0 + 5*sg;      % for plotting
a = a0 - 5*sg;
x = linspace(a, b, N+1);

%% Analogous Normal Distribution
  dnm = sg *sqrt(2*pi);       % denominator
  rt = (x - mu) /sg;         % ratio
f = 1/dnm * exp(-0.5 *rt.^2);

%% Analogous Cumulative Function
    rt_c = (x - mu) /(sg*sqrt(2));   % ratio
fc = 0.5* (1+ erf(rt_c) );

%% From Model
P_o1 = 7.15;        % outlet pressure; MPa
T_av1 = 284.9;      % Avg temp; dC
m_fl1 = 20.3;       % mass inflow; t/hr
wd = 132.5;         % assembly width; mm
rd = 8.0;           % assembly radius; mm
cf = 1.0;           % friction coeff
ksp = 1.2;          % spacer-grid p-loss factor

mu2 = bfbt_T(P_o1, T_av1, m_fl1, wd, rd, cf, ksp); 

    %% Plotting
% figure(1)
% plot(x,f), grid on, hold on
% plot(x,fc)
% xlabel('Values')
% ylabel('Probability Dist.')

%% For Sample
X = sort(DP_e);         % radom variable
mu_s = mean(DP_e);      % sample mean
sg_s = std(DP_e);       % sample STD
N = length(X);

  sg_z = sg_s/sqrt(N);    % STD-z
Z = (mu_s - mu2) /sg_z % new RV

    %% Plot



%% Functions ####################
function y = std_ud(X)
    N = length(X);
    mu = sum(X)/N;  % mean finder
       
    sqsm = 0;       % square sum
   for i=1:1:length(X)
        sqsm = sqsm + (X(i) - mu)^2;
    end 
    
    y = sqrt(sqsm/(N-1));
end

