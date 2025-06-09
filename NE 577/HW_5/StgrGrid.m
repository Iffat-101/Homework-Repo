function [x_ou, x_ev, y_ov, y_eu, x_ep, y_ep, dx, dy] = StgrGrid(nx, ny, Lx, Ly)
    % x_ou = x, odd, @ position of u; x_ev = x, even, @ position of v;
    % y_ov = y, odd, @ position of v; y_eu = y, even, @ position of u;
    % x_ep = x, even, @ position of p
    % x_ep = y, even, @ position of p

%  global nx ny Lx Ly
%% Grid
    % defining grids, to put values of x and y
    x0 = linspace(0, nx, 2*nx +1);  % additional 1 for boundaries
    x_ou = x0(1:2:end-1);       % excluding the right boundary
    x_ev = x0(2:2:end);         % takes up even values
    dx = x_ev(2) - x_ev(1);

    y0 = linspace(0, ny, 2*ny +1); % additional 1 for boundaries
    y_ov = y0(1:2:end-1);       % excluding the top boundary
    y_eu = y0(2:2:end);         % takes up even values
    dy = y_ov(2) - y_ov(1);
  
    x_ep = x_ev;    % defining grid for p's
    y_ep = y_eu;
end
