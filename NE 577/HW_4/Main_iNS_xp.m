% function Main_iNS
clear all, close all, clc
    global np nu nx ny dx dy dt Lx Ly iU iV iP...
        uBC_L uBC_R uBC_B uBC_T vBC_L vBC_R vBC_T vBC_B

    %% Constants, Grid
    nx = 100;    ny = 30;
    Lx = 0.04;   Ly = 0.02;
    np = nx*ny;
    nu = (nx-1)*ny + nx*(ny-1);
    
 %% Pointers
    [iP, iU, iV] = GenPointer(nx, ny);

 %% Grid
    % defining grids, to put values of x and y
      % This's no bbig deal for numerical solution. Only needed for analytical testing
    x0 = linspace(0, Lx, 2*nx +1); % additional 1 for boundaries
    x_ou = x0(1:2:end-1); % excluding the right boundary
    x_ev = x0(2:2:end); % takes up even values
    dx = x_ev(2) - x_ev(1);

    y0 = linspace(-Ly/2, Ly/2, 2*ny +1); % additional 1 for boundaries
    y_ov = y0(1:2:end-1); % excluding the top boundary
    y_eu = y0(2:2:end); % takes up even values
    dy = y_ov(2) - y_ov(1);
    
%     x_ep = x_ev; % defining grid for p's
%     y_ep = y_eu;

%% Boundary Conditions
% uBC_L = -185 *y_eu.^2 + 0.0185;
uBC_L = 0.0125  .*ones(size(y_eu));
uBC_R = uBC_L;  % Periodic BC
uBC_B = 0  .*ones(size(x_ou));
uBC_T = 0  .*ones(size(x_ou)); 

vBC_L = 0  .*ones(size(y_ov));
vBC_R = 0  .*ones(size(y_ov));
vBC_B = 0  .*ones(size(x_ev));
vBC_T = 0  .*ones(size(x_ev));

%% Initialize the Flow, Inhomogenous IC
    [P_now, U_now_i] = Flow_Initial(nx, ny);
        % P doesn't need initialization
    U_old_i = zeros(size(U_now_i));   % needed for advection
   
 %% Constants, Flow
    dt = 1e-5; % verify
    neu = 1e-6;
   
 %% time stepping using fractional step
 U_old = U_old_i;
 U_now = U_now_i;   % loop initials

 lim = 5e2; % 1e10; % iteration stopper, check
 tc = 1; % iteration counter
 for it = 1:1:lim
    %% Fractional Step: Stage 1
    % Defining operator S
    Lp_now = OpLap(U_now); % putting laplacian operator
    Su_now = U_now + dt/2 *neu *Lp_now;  % make it into matrix form
    
    % Defining advection vectors
    A_old = OpAdv(U_old, uBC_L, uBC_R, uBC_B, uBC_T, vBC_L, vBC_R, vBC_T, vBC_B);
    A_now = OpAdv(U_now, uBC_L, uBC_R, uBC_B, uBC_T, vBC_L, vBC_R, vBC_T, vBC_B);
    
    % Defining BC_Laplace
    bc_L_new = OpBC_Lap(uBC_L, uBC_R, uBC_B, uBC_T, vBC_L, vBC_R, vBC_T, vBC_B);
    
    % RHS Matrix
    RHS_b1 = Su_now + dt/2 *(3*A_now - A_old) + dt *neu *bc_L_new; % two are same
        
    st1 = 1e-6;     % default
    u_F = CGs_1(RHS_b1, nu, st1, dt, neu); % beaware of size
    % operator R is inside the CGS_1 solver
        
    %% Fractional Step: Stage 2
    Div_uF = OpDiv(u_F);
    bc_D_new = OpBC_Div(uBC_L, uBC_R, vBC_B, vBC_T);
    
    RHS_b2 = 1/dt *Div_uF + 1/dt *bc_D_new;
    
    st2 = 1e-6;     % default
    P_new = CGs_2(RHS_b2, np, st2); % beware size
      % DG operator is inside the CGs_2 solver with assumed: R^-1 = I
       
    %% Fractional Step: Stage 3
    Grd_P_new = OpGrad(P_new);
    U_new = u_F - dt *Grd_P_new;
        
       %% Values Reassignments
    U_old = U_now;
    U_now = U_new;
    tc = tc +1; % iteration counter
  
      %% Contour Movie
%         if mod(it, 10) == 0
%             chrlp = 'Contour, u'; chrlp2 = 'Contour, v';
%                 % printing the characters
%             PrintGenVct(U_new, 11, chrlp, chrlp2)
%         end
    
        %% Loop Breaker
          crit_up = norm(abs(U_now - U_old), inf); % breaking criterion
          crit_dn = dt; % norm(abs(U_now), inf);
        crit = crit_up/crit_dn;
        tol = 1e-9;
        if ( crit < tol ) %# Criteria, very Vital
            break;         % if more than 1e-6, loop stops quick 
        end
 end

    %% Matrix Viewer, at loop end
%         % this section is to see U_old and U_now assuming values at end of
%         % each iteration
%     [U_old_u, U_old_v] = VectToMat(U_old);
%     [U_old_u_vw, U_old_v_vw] = MatView(U_old_u, U_old_v);
%         % makes the Mat' look like how it's seen in diagram
%     [U_now_u, U_now_v] = VectToMat(U_now);
%     [U_now_u_vw, U_now_v_vw] = MatView(U_now_u, U_now_v);

    %% Appending BCs, and making Mat'
    % for u's
    X_sol = U_new ;
    [Mat_u, Mat_v] = BCAppend(X_sol);
        % creating 2 BC-included Mat' from solution vector
        
    % for P's
    [Mat_p] = VectToMatP(P_new);
    
    % for Grad-P
    [Mat_pG] = VectToMatP(Grd_P_new);

    %% Plotting, at x-y plane
    % for u's
%     chtit1 = 'Countour, $$ u_1 $$';
%     chrx1 = 'Domain length, $$ x_1 $$'; chry1 = 'Domain Width, $$ x_2 $$';
%     [Vw_u] = PrintMat(Mat_u, 11,  chtit1, chrx1, chry1);
    
    % for v's
%     chtit2 = 'Countour, $$ v_1 $$';
%     chrx2 = 'Domain length, $$ x_1 $$'; chry2 = 'Domain Width, $$ x_2 $$';  % for v's
%     [Vw_v] = PrintMat(Mat_v, 12,  chtit2, chrx2, chry2); 
    
    % for p's
    chtitP = 'Countour, $$ P^{n+1} $$';
    chrxP = 'Domain length, $$ x_1 $$'; chryP = 'Domain Width, $$ x_2 $$';
    [Vw_p] = PrintMat(Mat_p, 13,  chtitP, chrxP, chryP);
    
    % for Grd-p's
    chtitPG = 'Countour, $$ \nabla P^{n+1} $$';
    chrxPG = 'Domain length, $$ x_1 $$'; chryPG = 'Domain Width, $$ x_2 $$';
    [Vw_Grd_p] = PrintMat(Mat_pG, 14,  chtitPG, chrxPG, chryPG);

        %% Curve Plot
    [y_eu_anl, u1_y_anl] = PlotGenAnlytic(40);
        % Plots the u value Analytic Function

    [y_eu_nmr, u1_y_nmr] = PlotGenNmr(X_sol, 40);
        % Plots the u value from Numerical Results

% end
