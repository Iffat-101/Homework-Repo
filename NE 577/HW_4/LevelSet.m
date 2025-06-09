% function qo = LevelSet(qi)
clear all, close all, clc
    global np nu nx ny dx dy Lx Ly iU iV iP...
%         cx nx_rd ny_rd x_ep y_ep LS
    %% Constants
    nx = 1000;   ny = 500; % changeable 
    Lx = 10;    Ly = 5; % constant
%     np = nx*ny;
%     nu = (nx-1)*ny + nx*(ny-1);
    
    dia = Ly/10; rad = dia/2;
    n_bb = 5; n_dp = 5; % no. of bub/drplet in domain

    ny_rd = ny/Ly *rad; % resoltn of rad', along y
    ny_rd = round(ny_rd); % rounding up rad
    
    nx_rd = nx/Lx *rad; % resoltn of rad', along x 
    nx_rd = round(nx_rd);
    
 %% Pointers
%     [iP, iU, iV] = GenPointer(nx, ny);

 %% Grid
 [x_ou, x_ev, y_ov, y_eu, x_ep, y_ep, dx, dy] = StgrGrid(nx, ny, Lx, Ly);
    % instead of lenght, it takes resolutin

    %% Initialize output
    LS = nan(nx, ny);

  %% Defining Liquid and Vapor
    liq = 2;
    ny_lq = ny/Ly *2;   % resoltn', along liq' height
    ny_lq = round(ny_lq); % rounding up the value

    for i = 1:1:nx
        for j = 1:1: ny_lq
            LS(i,j) = 1;  % intakes value for liq'
        end
       
        for j = ny_lq+1 :1 :ny
            LS(i,j) = 0;  % intakes value for gas
        end
    end
    
%% Bubbles in Liquid
Lx_c1 = linspace(1, Lx-1, n_bb); % 5 number of bubbles
Ly_c1 = linspace(1, Ly-1, n_bb); % won't be needed

% cx1 = ny/Lx *[50, 150, 250, 350, 450];
cx1 = nx/Lx *Lx_c1;
    cx1 = round(cx1); % rounding up the value for lesser resolution
cy1 = [50, 150,  90, 175,  75];  %## Input

for k = 1:1:5
%     cx(k) = 50; cy(k) = 50;   % center indice
    if cy1(k) > ny_lq
        error('Bubble %d: value of cy has to be less than %d', k, ny_lq)
    end 
      
    for i = (cx1(k) -nx_rd):1: (cx1(k) +nx_rd)
        for j = (cy1(k) -nx_rd):1: (cy1(k) +ny_rd)
            val = ( sqrt( (x_ep(i) -cx1(k)).^2 + (y_ep(j) -cy1(k)).^2 ) - nx_rd ) /nx_rd;
    %         val = ( sqrt( (x_ep(i) -Lx_c1).^2 + (y_ep(j) -Ly_c1).^2 ) - nx_rd ) /nx_rd;
            if val < 0
                LS(i,j) = val;
            else
                LS(i,j) = 1;
            end
        end
    end 
end 
 
%% Droplets in Vapor
Lx_c2 = linspace(1, Lx-1, n_dp);
Ly_c2 = linspace(1, Ly-1, n_dp); % Won't be needed

cx2 = nx/Lx *Lx_c2;
    cx2 = round(cx2); % rounding up the value for lesser resolution
cy2 = [300, 280, 400, 475, 275]; %## Input

for k = 1:1:5
%     cx(k) = 50; cy(k) = 50;   % center indice
    if cy2(k) < ny_lq
        error('Droplet %d: value of cy2 has to be greater than %d', k, ny_lq)
    end 
  
    
    for i = (cx2(k) -nx_rd):1: (cx2(k) +nx_rd)
        for j = (cy2(k) -nx_rd):1: (cy2(k) +ny_rd)
            val = ( nx_rd - sqrt( (x_ep(i) -cx2(k)).^2 + (y_ep(j) -cy2(k)).^2 )) /nx_rd;
            if val > 0
                LS(i,j) = val;
            else
                LS(i,j) = 0;
            end
        end
    end 
end 
LS_vw = flipud(LS'); % for viewing level set as diag

%% Plotting
tit1 = 'Level Set';
chx1 = 'Domain Length, x'; chy1 = 'Domain Width, y';
PrintMat(LS, 11, tit1, chx1, chy1);
% axis([1 10 1 5]);

% end
