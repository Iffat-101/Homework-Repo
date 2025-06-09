clear all, close all, clc

w = linspace(0,1,100);
 l = length(w);

ph_r = 1;       % constant radial flux-density
C = ph_r/(4*pi);    % constant

ph_0 = C *ones(l);
ph_1 = C *w;
ph_2 = C *0.5 *(3*w.^2 - 1);
ph_3 = C *0.5 *(5*w.^3 - 3.*w);
ph_4 = C *0.125 *(35*w.^4 - 30.*w.^2 +3);

figure(1)
    grid on, hold on 
plot(w, ph_0)
plot(w, ph_1)
plot(w, ph_2)
plot(w, ph_3)
plot(w, ph_4)
    legend
legend('\phi_0', '\phi_1', '\phi_2', '\phi_3', '\phi_4')
    hold off

xlabel('\omega')
ylabel('\phi')
title('Flux Density (\phi) vs. Direction Cosine (\omega)')



