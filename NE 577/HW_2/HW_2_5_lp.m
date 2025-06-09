% clear all, close all, clc
%% #################### HW-2.5 ##################################
%% Section abc
N_c = 8730112;              % total no. of cores

%% No. of Processable Nodes (Prob d)
Cp = 3.25 *10^12;
Cp_Nd = Cp /10^6;
Cp_mx = 1.102 *10^18;
N_dd = Cp_mx /Cp_Nd;    % N_d for prob. d

%% Calculation
a_v = [0.25, 0.5, 0.75, 1];
n_o = [2048, 4096, 8192];

    %% No. of Proocessable Nodes
N_d = [];   % initialize N_d
N_a = N_c .*a_v;
len1 = length(a_v);
for i = 1:1:len1
    N_dL = N_a(i) .*n_o;          % for the loop
    N_d = [N_d, N_dL, N_dd];      % appending N_d after each loop-work
                      % here appends an extra N_dd after each loop
end

    %% Total volume of bubble
len2 = length(N_d);
v_f = [0.02, 0.04, 0.08];
V_bt = [];              % initializing V_bt
for j = 1: 1: len2      % need more smart
    V_btL = N_d(j) .*v_f;
    V_bt = [V_bt, V_btL];    % Appending after each loop
end

    %% No. of bubbles
len3 = length(V_bt);
n_d = [10, 20, 30];     % no. of bubble acrs dia
n_bb = [];
for k = 1:1:len3
    n_bbL = 6 *V_bt(k) ./(n_d.^3 *pi);
                        % this represents a matrix
    n_bb = [n_bb, n_bbL];
end

%% Plotting
    %% No. of processable Nodes
figure(1)
len_a = length(N_d);
x_a = linspace(1, len_a, len_a); % increase by 1
y_a = N_d;
plot(x_a,y_a, 'k'), grid on
xlabel('Serial'); ylabel('No of Processable Nodes')
title('Serial vs No. of processable nodes')

    %% Total Volume of bubble
figure(2)
len_b = length(V_bt);
x_b = linspace(1, len_b, len_b);
y_b = V_bt;
plot(x_b, y_b, 'g'); grid on
xlabel('Serial'); ylabel('Volume of bubble')

    %% No. of bubbles
figure(3)
len_c = length(n_bb);
x_c = linspace(1,len_c,len_c);
y_c = n_bb;
plot(x_c, y_c), grid on
xlabel('Serial'); ylabel('No. of bubbles')
title('Serial vs No. of bubbles')

