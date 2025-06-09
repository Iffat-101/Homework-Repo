function [X, U] = PlotGenAnlytic(fn)
% Plots the analytic function
%% Analytical Funtion
global nx ny dx dy Lx Ly

 %% Grid
    y0 = linspace(-Ly/2, Ly/2, 2*ny +1); % additional 1 for boundaries
%     y_ov = y0(1:2:end-1); % excluding the top boundary
    y_eu = y0(2:2:end); % takes up even values
    y_eu_ex = [-Ly/2, y_eu, Ly/2]; % appended boundary nodes
    
    u1_y = -185 *y_eu.^2 + 0.0185;
    u1_y_ex = [0, u1_y, 0]; % appended boundary node values
    
    %% The Plot
    figure(fn)
    plot(y_eu_ex, u1_y_ex, 'LineWidth', 1.5), grid on, hold on
    xlabel('Channel Width, $$ x_2(m) $$', 'Interpreter', 'Latex')
    ylabel('Velocity, $$ u_1(m/s) $$', 'Interpreter', 'Latex')
    title(' $$ x_2 $$ vs $$ u_1 $$', 'Interpreter', 'Latex')
    
    %% Output
    X = y_eu_ex; U = u1_y_ex; 
end
