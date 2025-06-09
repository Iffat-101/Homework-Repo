clear all, close all, clc
%% #################### HW-2.2 ##################################
%% Part a
T = 320;
D = 0.007; r = D/2;
kp = 2/r; % curvature 
gm = 0.07275 *(1 - 0.002 *(T- 291));

dl_P = gm * kp       ;% pressure drop
tht0 = 0: pi/180: 2*pi;

    %% Plotting
kp0 = kp *ones(size(tht0));
r0 = 2 ./kp0;

figure(1)
% ax = polaraxes;     % declares a var for polar-axes
polarplot(tht0, r0, 'g', 'linewidth', 2), hold on
rlim([0 0.004])

% ax.ThetaZeroLocation = 'left';  % puts tht position at left
title('Bubble shape at 320 K')
 
%% Part c-1
tht1 = 0: pi/180: pi;
tht2 = pi: pi/180: 2*pi;    % declaring the range
        % start @0, increse 1 dg, upto 180 dg

T_1 = 270; T_2 = 320;
T_vr1 = T_1 + (T_2 - T_1) .*tht1/pi;   % declares 1st range 
T_vr2 = T_2 - (T_2 - T_1) .*(tht2 - pi)/pi; % declare 2nd range

gm1 = 0.07275 *(1 - 0.002 .*(T_vr1 - 291));
gm2 = 0.07275 *(1 - 0.002 .*(T_vr2 - 291));

kp_b1 = dl_P ./gm1;
kp_b2 = dl_P ./gm2;

    %% Plotting
tht = [tht1, tht2];
kp_b = [kp_b1, kp_b2];
r_b = 2 ./kp_b;

figure(2)
polarplot(tht, r_b, 'b', 'linewidth', 1)
rlim([0 0.004])
title('Bubble shape at 270 K - 320 K')

%% Part c-2
tht1_2 = 0: pi/180: pi;
tht2_2 = pi: pi/180: 2*pi;    % declaring the range
        % start @0, increse 1 dg, upto 180 dg

T_1_2 = 270; T_2_2 = 370; % the change is here 
T_vr1_2 = T_1_2 + (T_2_2 - T_1_2) .*tht1_2/pi;
T_vr2_2 = T_2_2 - (T_2_2 - T_1_2) .*(tht2_2 - pi)/pi;

gm1_2 = 0.07275 *(1 - 0.002 .*(T_vr1_2 - 291));
gm2_2 = 0.07275 *(1 - 0.002 .*(T_vr2_2 - 291));

kp_b1_2 = dl_P ./gm1_2;
kp_b2_2 = dl_P ./gm2_2;

    %% Plotting
tht_2 = [tht1_2, tht2_2];   % Appending
kp_b_2 = [kp_b1_2, kp_b2_2];
r_b_2 = 2 ./kp_b_2;

figure(3)
polarplot(tht_2, r_b_2, 'r', 'linewidth', 1)
rlim([0 0.004])
title('Bubble shape at 270 K - 370 K')

