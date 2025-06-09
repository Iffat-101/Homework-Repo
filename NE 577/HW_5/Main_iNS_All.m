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
   
        %% Loop Breaker
          crit_up = norm(abs(U_now - U_old), inf); % breaking criterion
          crit_dn = dt; % norm(abs(U_now), inf);
        crit = crit_up/crit_dn;
        tol = 1e-9;
        if ( crit < tol ) %# Criteria, very Vital
            break;         % if more than 1e-6, loop stops quick 
        end
 end

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

%% ########## Rest of the functions are plotted Alphabetically
%% ##########   ############################################################
function [X, Y] = BCAppend(Vec)
    % Vec = input Vector; nx, ny = matrix size; fn = figure no.; chr input
    % first, clip the vector into two; generate matrix > contour for 'em
    
    global uBC_L uBC_R uBC_B uBC_T vBC_L vBC_R vBC_T vBC_B
    [Mat_u, Mat_v] = VectToMat(Vec);
    
    Mat_u(1,:) = uBC_L; % Replacing the BCs
    Mat_v(:,1) = vBC_B;

    Mat_uT = Mat_u';  % transpose the Matrix makin it similar to diagrm
    Mat_uT = [Mat_uT, uBC_R']; % appending right BC
    uBC_Te = [uBC_T, 0]; % adding extra 0 to matchup size
    uBC_Be = [uBC_B, 0]; 
    Mat_uT = [uBC_Be; Mat_uT; uBC_Te]; % Appending the boundaries
        % these are half-grid appendix, so, makes domain look longer
    
    Mat_vT = Mat_v';
    Mat_vT = [Mat_vT; vBC_T];
    vBC_Le = [vBC_L, 0];
    vBC_Re = [vBC_R, 0];
    Mat_vT = [vBC_Le', Mat_vT, vBC_Re'];

    X = Mat_uT'; % reverting back to Matlab form
    Y = Mat_vT';

    Vw_1 = flipud(Mat_uT); % for viewing how it looks like in mat form
    Vw_2 = flipud(Mat_vT);  
end

%% ##########	############################################################
function X_f = CGs_1(RHS_b, sz, sc,  dt, neu)
    % RHS_b = RHS vector in the equation
    % sz = grid-size of the domain (nu/np)
    % sc = stopping criterion
%% Initialization
X = zeros(sz,1);
P = zeros(sz,1);
R = RHS_b;

%% Main Loop
for k1=1:1:1000
    k = k1 - 1; % so, starts with k = 0;
    if (k == 0) % valid only for 1st iteration
        bet = 0;
    else      % row * col' = num; col' * row = matrix
        bet = (R'*R)/(R_ld' *R_ld); % bet, a number
    end 
    
    P = R + bet * P;
    Ap =  OpLap(P);
    %%# Defining Operator R
    Ap = Ap - 1/2 *dt *neu *Ap;
    
    alf = (R'*R)/(P' *(Ap)); % a constant
    X = X + alf * P;        % solution vct
    R_ld = R;    % previous value assigned to R_ld
    Ax =  OpLap(X);
    Ax = Ax - 1/2 * dt *neu *Ax;
    R = RHS_b - Ax; % new value assigned to R
    
    if (norm(R)/norm(RHS_b) <= sc)
        break;
    end
    
    if (k == (sz+1)) % Check, extra
        error('CG not working')
    end
end
%
X_f = X;

end

%% ##########	############################################################
function X_f = CGs_2(RHS_b, sz, sc)
    % RHS_b = RHS vector in the equation
    % sz = grid-size of the domain (nu/np)
    % sc = stoppiing criterion
%% Initialization
X = zeros(sz,1);
P = zeros(sz,1);
R = RHS_b;

%% Main Loop
for k1=1:1:1000
    k = k1 - 1; % so, starts with k = 0;
    if (k == 0) % valid only for 1st iteration
        bet = 0;
    else
        bet = (R'*R)/(R_ld' *R_ld);
    end
    
    P = R + bet * P;
    %%# Applying the 'DG' operator
    Ap = OpGrad(P);
    Ap = OpDiv(Ap);

    alf = (R'*R)/(P' *(Ap)); 
    X = X + alf * P; 
    R_ld = R;    % assigning prev R
    %%# Applying the DG operator
    Ax = OpGrad(X);
    Ax = OpDiv(Ax);
    R = RHS_b - Ax; % new R

    if (norm(R)/norm(RHS_b) <= sc)
        break;
    end
    
    if (k == 9999)  % Check, extra added
        error('CG not working')
    end
end
%
X_f = X;

end

%% ##########	############################################################
function [P_now, U_now_i] = Flow_Initial(nx, ny)
    % takes the grid-size as input and assigns initial values for P & U
    % returns long vector (u+v) for U and a vector for P
    
%   nx = 5; ny = 5;
   %% Memory Allocation
   iP = nan(nx, ny);
   iU = nan(nx, ny);
   iV = nan(nx, ny);
   
      %% initialize values
    Uu_o1_Mt = iU; Uv_o1_Mt = iV; P_o1_Mt = iP;
    
   %% for U's
       % X(u) component
   for i = 2:1:nx
        for j = 1:1:ny
            Uu_o1_Mt(i,j) = 0; % no pre-u in the domain
        end 
    end
    
        % Y(v) component
    for i = 1:1:nx
        for j = 2:1:ny
            Uv_o1_Mt(i,j) = 0; % no pre-v in the domain
        end 
    end
    
    % Convering into Long Vector
    [Uu_o1L, Uv_o1L] = MatToVect(Uu_o1_Mt, Uv_o1_Mt);
    U_now_i = [Uu_o1L, Uv_o1L]';   % Appending two vectors, making col'
        % Because operators deal with u+v long vectors

    %% for P's
    for i = 1:1:nx
        for j = 1:1:ny
            P_o1_Mt(i,j) = 0; % no pre-P in the cavity
        end 
    end
    
    % Converting into Long Vector
    P_now = MatToVectP(P_o1_Mt);
    P_now = P_now'; % making col'

end

%% ##########	############################################################
function [iP, iU, iV] = GenPointer(nx, ny)
% global nx ny
%% Memory Allocation
    iP =  nan(nx, ny);
    iU =  nan(nx, ny); 
    iV =  nan(nx, ny);  % allocating NaN for variables
    
    %% Pointer matrix for P
    id_p = 1; % index to be used in vector variable P
    for i =1:1:nx
        for j = 1:1:ny
            iP(i,j) = id_p; % putting number in pointers
            id_p = id_p + 1;       
        end
    end 
%     iP; % output

    %% Pointer Matrix for ux
    id_u = 1; % index to be used in vector variable u = [ux; uy]
    for i =2:1:nx   % starts from 2
        for j = 1:1:ny
            iU(i,j) = id_u;
            id_u = id_u + 1;
        end
    end
%     iU; % output

    %% Pointer Matrix for vy
    for i =1:1:nx
        for j = 2:1:ny  % starts from 2
            iV(i,j) = id_u;
            id_u = id_u + 1;   % long vector index    
        end
    end
%     iV; % output

    %% Visualizer
    iP_vw = flipud(iP'); % viewing the memory locations as like diagram
    iU_vw = flipud(iU');
    iV_vw = flipud(iV');
end


%% ##########	############################################################
function [Vec_u, Vec_v] = MatToVect(Mat_u, Mat_v)
    global nx ny np nu 
        % here, boundary values are clipped off for u and v
        % then Matrix is converted into row vector
    %% Memory Allocation
    V_uL = nan(1, (nx-1)*ny );
    V_vL = nan(1, nx *(ny-1));

    %% Creating Long vector
    id_u = 1;
    for i = 2:1:nx
        for j = 1:1:ny
            V_uL(id_u) = Mat_u(i,j);
            id_u = id_u +1;
        end % generating long-vector, excluding the boundary values
    end
    Vec_u = V_uL;

     % craeting long vector with for v
     id_v = 1;
    for i = 1:1:nx
        for j = 2:1:ny
            V_vL(id_v) = Mat_v(i,j);
            id_v = id_v +1;
        end 
    end
    Vec_v = V_vL;
end 


%% ##########	############################################################
function Vec_p = MatToVectP(Mat_p)
    global nx ny nu np 
        % here, boundary values are 'Not' clipped off for u and v
        % then Matrix is converted into row vector
    %% Memory Allocation
    V_pL = nan(1, np);

    %% Creating Long vector
    id_p = 1;
    for i = 1:1:nx  % retains 1st node values
        for j = 1:1:ny
            V_pL(id_p) = Mat_p(i,j);
            id_p = id_p +1;
        end
    end
    Vec_p = V_pL;

end


%% ##########	############################################################
function qo = OpAdv (qi, uBC_L, uBC_R, uBC_B, uBC_T, vBC_L, vBC_R, vBC_T, vBC_B)
    global nu iU iV nx ny dx dy 
    % input : u type n v type (mixed)
    % output: u type
    
    %% Input size check
    if size (qi, 2) ~= 1
        error ('Adv:: Need a colun vector input!!')
    end 
    
    if size (qi, 1) ~= nu
        error ('Adv:: input vector size mismatch!!')
    end
    
    %% Initialize output
    qo = nan(nu, 1);
    
 %% 1. u-component
    %% Inner Domain
    for i = 3:1:(nx-1)
        for j = 2:1:(ny-1)     % defining u_br_x                        /squaring it
            qo(iU(i, j)) = - ( - ( qi(iU(i-1, j)) + qi(iU(i, j))   ) /2 * ( qi(iU(i-1, j))   + qi(iU(i, j))   ) /2       ... % u_x_lf
                               + ( qi(iU(i, j))   + qi(iU(i+1, j)) ) /2 * ( qi(iU(i, j))     + qi(iU(i+1, j)) ) /2 ) /dx ... % u_x_rt                             
                           - ( - ( qi(iU(i, j-1)) + qi(iU(i, j))   ) /2 * ( qi(iV(i-1, j))   + qi(iV(i, j))   ) /2       ... % bt
                            ...     % u_y_bt                             v_x_bt
                               + ( qi(iU(i, j))   + qi(iU(i, j+1)) ) /2 * ( qi(iV(i-1, j+1)) + qi(iV(i, j+1)) ) /2 ) /dy ;   % tp
                                    % u_y_tp                            v_x_tp            
        end
    end
    
    %% Edges
    % left inner (T1)
    i = 2;
        for j = 2:1:(ny-1)         % only change here
            qo(iU(i, j)) = - ( - ( uBC_L(j)       + qi(iU(i, j))   ) /2 * ( uBC_L(j)         + qi(iU(i, j))   ) /2       ...
                               + ( qi(iU(i, j))   + qi(iU(i+1, j)) ) /2 * ( qi(iU(i, j))     + qi(iU(i+1, j)) ) /2 ) /dx ...                               
                           - ( - ( qi(iU(i, j-1)) + qi(iU(i, j))   ) /2 * ( qi(iV(i-1, j))   + qi(iV(i, j))   ) /2       ... 
                               + ( qi(iU(i, j))   + qi(iU(i, j+1)) ) /2 * ( qi(iV(i-1, j+1)) + qi(iV(i, j+1)) ) /2 ) /dy ;            
        end
    
    % right inner (T1)
    i = nx;
        for j = 2:1:(ny-1)
            qo(iU(i, j)) = - ( - ( qi(iU(i-1, j)) + qi(iU(i, j)) )   /2 * ( qi(iU(i-1, j))   + qi(iU(i, j))   ) /2       ...
                               + ( qi(iU(i, j))   + uBC_R(j)     )   /2 * ( qi(iU(i, j))     + uBC_R(j)       ) /2 ) /dx ...                            
                           - ( - ( qi(iU(i, j-1)) + qi(iU(i, j)) )   /2 * ( qi(iV(i-1, j))   + qi(iV(i, j))   ) /2       ... 
                               + ( qi(iU(i, j))   + qi(iU(i, j+1)) ) /2 * ( qi(iV(i-1, j+1)) + qi(iV(i, j+1)) ) /2 ) /dy ;    
        end
     
    % bottom innner (T2)
    j = 1;
         for i = 3:1:(nx-1)
            qo(iU(i, j)) = - ( - ( qi(iU(i-1, j)) + qi(iU(i, j))   ) /2 * ( qi(iU(i-1, j))   + qi(iU(i, j))   ) /2       ...
                               + ( qi(iU(i, j))   + qi(iU(i+1, j)) ) /2 * ( qi(iU(i, j))     + qi(iU(i+1, j)) ) /2 ) /dx ...                            
                           - ( - ( uBC_B(i)       + uBC_B(i)       ) /2 * ( vBC_B(i-1)       + vBC_B(i)       ) /2       ... % note vBC
                               + ( qi(iU(i, j))   + qi(iU(i, j+1)) ) /2 * ( qi(iV(i-1, j+1)) + qi(iV(i, j+1)) ) /2 ) /dy ;    
         end

    % top inner (T2)
    j = ny;
        for i = 3:1:(nx-1)
            qo(iU(i, j)) = - ( - ( qi(iU(i-1, j)) + qi(iU(i, j))   ) /2 * ( qi(iU(i-1, j))   + qi(iU(i, j))   ) /2       ...
                               + ( qi(iU(i, j))   + qi(iU(i+1, j)) ) /2 * ( qi(iU(i, j))     + qi(iU(i+1, j)) ) /2 ) /dx ...                            
                           - ( - ( qi(iU(i, j-1)) + qi(iU(i, j))   ) /2 * ( qi(iV(i-1, j))   + qi(iV(i, j))   ) /2       ... 
                               + ( uBC_T(i)       + uBC_T(i)       ) /2 * ( vBC_T(i-1)       + vBC_T(i)       ) /2 ) /dy ;    
        end
     
    %% Corners
    % left bottom
    i = 2;
    j = 1;
            qo(iU(i, j)) = - ( - ( uBC_L(j)       + qi(iU(i, j))   ) /2 * ( uBC_L(j)         + qi(iU(i, j))   ) /2       ...
                               + ( qi(iU(i, j))   + qi(iU(i+1, j)) ) /2 * ( qi(iU(i, j))     + qi(iU(i+1, j)) ) /2 ) /dx ...                            
                           - ( - ( uBC_B(i)       + uBC_B(i)       ) /2 * ( vBC_B(i-1)       + vBC_B(i)       ) /2       ... 
                               + ( qi(iU(i, j))   + qi(iU(i, j+1)) ) /2 * ( qi(iV(i-1, j+1)) + qi(iV(i, j+1)) ) /2 ) /dy ; 
    
    % right bottom
    i = nx;
    j = 1;
            qo(iU(i, j)) = - ( - ( qi(iU(i-1, j)) + qi(iU(i, j))   ) /2 * ( qi(iU(i-1, j))   + qi(iU(i, j))   ) /2       ...
                               + ( qi(iU(i, j))   + uBC_R(j)       ) /2 * ( qi(iU(i, j))     + uBC_R(j)       ) /2 ) /dx ...                            
                           - ( - ( uBC_B(i)       + uBC_B(i)       ) /2 * ( vBC_B(i-1)       + vBC_B(i)       ) /2       ... 
                               + ( qi(iU(i, j))   + qi(iU(i, j+1)) ) /2 * ( qi(iV(i-1, j+1)) + qi(iV(i, j+1)) ) /2 ) /dy ; 
    
 
    % left top
    i = 2;
    j = ny;
            qo(iU(i, j)) = - ( - ( uBC_L(j)       + qi(iU(i, j))   ) /2 * ( uBC_L(j)         + qi(iU(i, j))   ) /2       ...
                               + ( qi(iU(i, j))   + qi(iU(i+1, j)) ) /2 * ( qi(iU(i, j))     + qi(iU(i+1, j)) ) /2 ) /dx ...
                           - ( - ( qi(iU(i, j-1)) + qi(iU(i, j))   ) /2 * ( qi(iV(i-1, j))   + qi(iV(i, j))   ) /2       ... 
                               + ( uBC_T(i)       + uBC_T(i)       ) /2 * ( vBC_T(i-1)       + vBC_T(i)       ) /2 ) /dy ; 

    % right top
    i = nx;
    j = ny;
           qo(iU(i, j)) = - ( - ( qi(iU(i-1, j)) + qi(iU(i, j))   ) /2 * ( qi(iU(i-1, j))   + qi(iU(i, j))    ) /2       ...
                              + ( qi(iU(i, j))   + uBC_R(j)       ) /2 * ( qi(iU(i, j))     + uBC_R(j)        ) /2 ) /dx ...
                          - ( - ( qi(iU(i, j-1)) + qi(iU(i, j))   ) /2 * ( qi(iV(i-1, j))   + qi(iV(i, j))    ) /2       ... 
                              + ( uBC_T(i)       + uBC_T(i)       ) /2 * ( vBC_T(i-1)       + vBC_T(i)        ) /2 ) /dy ; 

 %% v-components
    %% Inner domain
      for i = 2:1:(nx-1)
        for j = 3:1:(ny-1)         % v_x_lf                                 % u_y_lf
            qo(iV(i, j)) = - ( - ( qi(iV(i-1, j)) + qi(iV(i, j))   ) /2 * ( qi(iU(i, j-1))   + qi(iU(i, j))   ) /2       ... % lf
                               + ( qi(iV(i, j))   + qi(iV(i+1, j)) ) /2 * ( qi(iU(i+1, j-1)) + qi(iU(i+1, j)) ) /2 ) /dx ... % rt
                               ...
                           - ( - ( qi(iV(i, j-1)) + qi(iV(i, j))   ) /2 * ( qi(iV(i, j-1))   + qi(iV(i, j))   ) /2       ... % bt
                               + ( qi(iV(i, j))   + qi(iV(i, j+1)) ) /2 * ( qi(iV(i, j))     + qi(iV(i, j+1)) ) /2 ) /dy ;   % tp
        end
      end
    
    %% Edges
    % left inner (T2)
     i = 1;
        for j = 3:1:(ny-1)
            qo(iV(i, j)) = - ( - ( vBC_L(j)       + vBC_L(j)       ) /2 * ( uBC_L(j-1)       + uBC_L(j)        ) /2       ...
                               + ( qi(iV(i, j))   + qi(iV(i+1, j)) ) /2 * ( qi(iU(i+1, j-1)) + qi(iU(i+1, j))  ) /2 ) /dx ...
                           - ( - ( qi(iV(i, j-1)) + qi(iV(i, j))   ) /2 * ( qi(iV(i, j-1))   + qi(iV(i, j))    ) /2       ...
                               + ( qi(iV(i, j))   + qi(iV(i, j+1)) ) /2 * ( qi(iV(i, j))     + qi(iV(i, j+1))  ) /2 ) /dy ;
        end
  
    % right inner  (T2)
    i = nx;
         for j = 3:1:(ny-1) 
            qo(iV(i, j)) = - ( - ( qi(iV(i-1, j)) + qi(iV(i, j))   ) /2 * ( qi(iU(i, j-1))  + qi(iU(i, j))     ) /2       ...
                               + ( vBC_R(j)       + vBC_R(j)       ) /2 * ( uBC_R(j-1)      + uBC_R(j)         ) /2 ) /dx ...
                           - ( - ( qi(iV(i, j-1)) + qi(iV(i, j))   ) /2 * ( qi(iV(i, j-1))  + qi(iV(i, j))     ) /2       ...
                               + ( qi(iV(i, j))   + qi(iV(i, j+1)) ) /2 * ( qi(iV(i, j))    + qi(iV(i, j+1))   ) /2 ) /dy ;
         end    

    % bottom inner (T1)
     j = 2; 
        for  i = 2:1:(nx-1)
            qo(iV(i, j)) = - ( - ( qi(iV(i-1, j)) + qi(iV(i, j))   ) /2 * ( qi(iU(i, j-1))   + qi(iU(i, j))    ) /2       ...
                               + ( qi(iV(i, j))   + qi(iV(i+1, j)) ) /2 * ( qi(iU(i+1, j-1)) + qi(iU(i+1, j))  ) /2 ) /dx ...
                           - ( - ( vBC_B(i)       + qi(iV(i, j))   ) /2 * ( vBC_B(i)         + qi(iV(i, j))    ) /2       ...
                               + ( qi(iV(i, j))   + qi(iV(i, j+1)) ) /2 * ( qi(iV(i, j))     + qi(iV(i, j+1))  ) /2 ) /dy ;
        end
    
     % top inner (T1)
     j = ny;
         for i = 2:1:(nx-1)
            qo(iV(i, j)) = - ( - ( qi(iV(i-1, j)) + qi(iV(i, j))   ) /2 * ( qi(iU(i, j-1))   + qi(iU(i, j))   ) /2       ...
                               + ( qi(iV(i, j))   + qi(iV(i+1, j)) ) /2 * ( qi(iU(i+1, j-1)) + qi(iU(i+1, j)) ) /2 ) /dx ... 
                           - ( - ( qi(iV(i, j-1)) + qi(iV(i, j))   ) /2 * ( qi(iV(i, j-1))   + qi(iV(i, j))   ) /2       ...
                               + ( qi(iV(i, j))   + vBC_T(i)       ) /2 * ( qi(iV(i, j))     + vBC_T(i)       ) /2 ) /dy ;
        end    

    %% Corners
    % left bottom
    i = 1;
    j = 2;
            qo(iV(i, j)) = - ( - ( vBC_L(j)       + vBC_L(j)       ) /2 * ( uBC_L(j-1)       + uBC_L(j)       ) /2       ...
                               + ( qi(iV(i, j))   + qi(iV(i+1, j)) ) /2 * ( qi(iU(i+1, j-1)) + qi(iU(i+1, j)) ) /2 ) /dx ... 
                           - ( - ( vBC_B(i)       + qi(iV(i, j))   ) /2 * ( vBC_B(i)         + qi(iV(i, j))   ) /2       ...
                               + ( qi(iV(i, j))   + qi(iV(i, j+1)) ) /2 * ( qi(iV(i, j))     + qi(iV(i, j+1)) ) /2 ) /dy ; 

    % right bottom
    i = nx;
    j = 2;
            qo(iV(i, j)) = - ( - ( qi(iV(i-1, j)) + qi(iV(i, j))   ) /2 * ( qi(iU(i, j-1))   + qi(iU(i, j))   ) /2       ...
                               + ( vBC_R(j)       + vBC_R(j)       ) /2 * ( uBC_R(j-1)       + uBC_R(j)       ) /2 ) /dx ... 
                           - ( - ( vBC_B(i)       + qi(iV(i, j))   ) /2 * ( vBC_B(i)         + qi(iV(i, j))   ) /2       ...
                               + ( qi(iV(i, j))   + qi(iV(i, j+1)) ) /2 * ( qi(iV(i, j))     + qi(iV(i, j+1)) ) /2 ) /dy ;
    
    % left top
    i = 1;
    j = ny;
            qo(iV(i, j)) = - ( - ( vBC_L(j)       + vBC_L(j)       ) /2 * ( uBC_L(j-1)       + uBC_L(j)       ) /2       ...
                               + ( qi(iV(i, j))   + qi(iV(i+1, j)) ) /2 * ( qi(iU(i+1, j-1)) + qi(iU(i+1, j)) ) /2 ) /dx ... 
                           - ( - ( qi(iV(i, j-1)) + qi(iV(i, j))   ) /2 * ( qi(iV(i, j-1))   + qi(iV(i, j))   ) /2       ...
                               + ( qi(iV(i, j))   + vBC_T(i)       ) /2 * ( qi(iV(i, j))     + vBC_T(i)       ) /2 ) /dy ;    

    % right top
    i = nx;
    j = ny;
            qo(iV(i, j)) = - ( - ( qi(iV(i-1, j)) + qi(iV(i, j))   ) /2 * ( qi(iU(i, j-1))   + qi(iU(i, j))   ) /2       ...
                               + ( vBC_R(j)       + vBC_R(j)       ) /2 * ( uBC_R(j-1)       + uBC_R(j)       ) /2 ) /dx ... 
                           - ( - ( qi(iV(i, j-1)) + qi(iV(i, j))   ) /2 * ( qi(iV(i, j-1))   + qi(iV(i, j))   ) /2       ...
                               + ( qi(iV(i, j))   + vBC_T(i)       ) /2 * ( qi(iV(i, j))     + vBC_T(i)       ) /2 ) /dy ;    
    
end


%% ##########	############################################################
function bcD = OpBC_Div(uBC_L, uBC_R, vBC_B, vBC_T)
    global np iP nx ny dx dy 
    %% BC Vector due to DiVergence
        % input:  BCs
        % output: p type (np elements)
    %% Initialize output
    bcD = zeros(np, 1);

    %% Edges
    % bottom inner
     j = 1;
         for i = 2:1:(nx-1)
             bcD(iP(i,j)) = (              0             ) /dx ... % instead of follows
                          + ( - vBC_B(i)                 ) /dy; % - qi(iV(i, j))
         end

    % top inner
    j = ny;
         for i = 2:1:(nx-1)
             bcD(iP(i,j)) = (              0             ) /dx ...
                          + (              + vBC_T(i)    ) /dy; % qi(iV(i, j+1))
         end

    % left inner
    i = 1;
         for j = 2:1:(ny-1)
             bcD(iP(i,j)) = ( - uBC_L(j)                ) /dx ... % - qi(iU(i, j))
                          + (              0            ) /dy;
         end 

    % right inner
    i = nx;
         for j = 2:1:(ny-1)
             bcD(iP(i,j)) = (               + uBC_R(j)   ) /dx ...  % qi(iU(i+1, j))
                          + (               0            ) /dy;
         end
    
    %% Corners
    % bottom left (need pinning)
    i = 1;
    j = 1;    % prof. gave it a hard-coded zero here; why? 
             bcD(iP(i,j)) = ( - uBC_L(j)                 ) /dx... % - qi(iU(i, j))
                          + ( - vBC_B(i)                 ) /dy;  % - qi(iV(i, j))
    % bottom right
    i = nx;
    j = 1;
             bcD(iP(i,j)) = (               + uBC_R(j)   ) /dx... % + qi(iU(i+1, j))
                          + (  - vBC_B(i)                ) /dy;  % - qi(iV(i, j))

    % top left
    i = 1;
    j = ny;
             bcD(iP(i,j)) = ( - uBC_L(j)                 ) /dx... % - qi(iU(i, j))
                          + (               + vBC_T(i)   ) /dy;  % + qi(iV(i, j+1))
    % top right
    i = nx;
    j = ny;
             bcD(iP(i,j)) = (               + uBC_R(j)    ) /dx... % + qi(iU(i+1, j))
                          + (               + vBC_T(i)    ) /dy;  % + qi(iV(i, j+1))


end 


%% ##########	############################################################
function bcL = OpBC_Lap(uBC_L, uBC_R, uBC_B, uBC_T, vBC_L, vBC_R, vBC_T, vBC_B)
    global nu iU iV  nx ny dx dy  
    % input: BCs
    % output: u-type (nu elements)
    
    %% input size check:
    % input BC's are all scalars
    
    %% Initialize output
    bcL = zeros (nu, 1);

  %% 1. u-Component
    %% Inner Domain
    % nothing here
    %% Edges
    % left inner    (type 1)
    i = 2;
        for j = 2:1:(ny-1)
            bcL(iU(i, j)) = ( + uBC_L(j)                                         ) /dx^2 ... 
                          + (                        0                           ) /dy^2 ;
        end

    % right inner   (type 1)
    i = nx;
        for j = 2:1:(ny-1)
            bcL(iU(i, j)) = (                                       + uBC_R(j)   ) /dx^2 ... 
                          + (                        0                           ) /dy^2 ;
        end

    % bottom inner  (type 2)
    j = 1;
        for i = 3:1:(nx-1) 
            bcL(iU(i, j)) = (                        0                           ) /dx^2 ...
                          + ( + 2*uBC_B(i)                                       ) /dy^2 ;
        end

    % top inner    (type 2)
    j = ny;
        for i = 3:1:(nx-1)
             bcL(iU(i, j)) = (                       0                           ) /dx^2 ...
                           + (                                      + 2*uBC_T(i) ) /dy^2 ; % 
        end 

    %% Corners
    % left bottom
    i = 2;
    j = 1;
            bcL(iU(i, j)) = ( + uBC_L(j)                                        ) /dx^2 ... 
                          + ( + 2*uBC_B(i)                                      ) /dy^2 ;
     
    % right bottom
    i = nx;
    j = 1;
            bcL(iU(i, j)) = (                                        + uBC_R(j) ) /dx^2 ... 
                          + ( + 2*uBC_B(i)                                      ) /dy^2 ;   

    % left top
    i = 2;
    j = ny;
            bcL(iU(i, j)) = ( + uBC_L(j)                                         ) /dx^2 ... 
                          + (                                       + 2*uBC_T(i) ) /dy^2 ; 

    % right top
    i = nx;
    j = ny;
            bcL(iU(i, j)) = (                                       + uBC_R(j)   ) /dx^2 ... 
                          + (                                       + 2*uBC_T(i) ) /dy^2 ; 

  %% 2. v-Component
    %% Inner Domain
    % nothing here
    %% Edges
    % left inner    (type 2)
    i = 1;
        for j = 3:1:(ny-1)
            bcL(iV(i, j)) = ( + 2*vBC_L(j)                                      ) /dx^2 ...
                          + (                          0                        ) /dy^2 ;
        end
      
    % right inner  (type 2) 
    i = nx;
        for j = 3:1:(ny-1)
            bcL(iV(i, j)) = (                                       + 2*vBC_R(j) ) /dx^2 ... 
                          + (                          0                         ) /dy^2 ;
        end
 
    % bottom inner  (type 1)
    j = 2;
        for i = 2:1:(nx-1)
            bcL(iV(i, j)) = (                          0                        ) /dx^2 ...
                          + ( + vBC_B(i)                                        ) /dy^2 ;
        end
 
    % top inner     (type 1)
    j = ny;
        for i = 2:1:(nx-1)
            bcL(iV(i, j)) = (                          0                        ) /dx^2 ...
                          + (                                       + vBC_T(i)  ) /dy^2 ;
        end

    %% Corners
    % left bottom
    i = 1;
    j = 2;
            bcL(iV(i, j)) = ( + 2*vBC_L(j)                                       ) /dx^2 ...
                          + ( + vBC_B(i)                                         ) /dy^2 ;
      
    % right bottom
    i = nx;
    j = 2;
            bcL(iV(i, j)) = (                                       + 2*vBC_R(j) ) /dx^2 ...
                          + ( + vBC_B(i)                                         ) /dy^2 ;

    % left top
    i = 1;
    j = ny;
            bcL(iV(i, j)) = ( + 2*vBC_L(j)                                       ) /dx^2 ...
                          + (                                       + vBC_T(i)   ) /dy^2 ;

    % right top
    i = nx;
    j = ny;
            bcL(iV(i, j)) = (                                       + 2*vBC_R(j)  ) /dx^2 ...
                          + (                                       + vBC_T(i)    ) /dy^2 ;

end

%% ##########	############################################################
function qo = OpDiv(qi)
    global np nu nx ny dx dy iU iV iP 
    %% Divergence operator
        % input:  u type (nu elements)
        % output: p type (np elements)
    
    %% Input size check
    if size(qi, 2) ~= 1
        error('Div:: input is not a column vector')
    end 

    if size(qi, 1) ~= nu
        error ('Div:: input is not a size mismatch')
    end 

    %% Initialize output
    qo = nan(np, 1);

    %% Inner Domain
     for i = 2:1:(nx-1)  % we exclude 2 edges here
         for j = 2:1:(ny-1)
             qo(iP(i,j)) = ( - qi(iU(i, j)) + qi(iU(i+1, j)) ) /dx...
                         + ( - qi(iV(i, j)) + qi(iV(i, j+1)) ) /dy;
         end 
     end
     
     %% Edges
     % bottom inner
     j = 1;
         for i = 2:1:(nx-1)
             qo(iP(i,j)) = ( - qi(iU(i, j)) + qi(iU(i+1, j)) ) /dx...
                         + (                + qi(iV(i, j+1)) ) /dy; % - qi(iV(i, j))
         end

     % top inner
     j = ny;
         for i = 2:1:(nx-1)
             qo(iP(i,j)) = ( - qi(iU(i, j)) + qi(iU(i+1, j)) ) /dx...
                         + ( - qi(iV(i, j))                  ) /dy; % qi(iV(i, j+1))
         end 

     % left inner
     i = 1;
         for j = 2:1:(ny-1)
             qo(iP(i,j)) = (                + qi(iU(i+1, j)) ) /dx... % - qi(iU(i, j))
                         + ( - qi(iV(i, j)) + qi(iV(i, j+1)) ) /dy;
         end    % the value in the blank goes to the BC vector

     % right inner
     i = nx;
         for j = 2:1:(ny-1)
             qo(iP(i,j)) = ( - qi(iU(i, j))                  ) /dx... % qi(iU(i+1, j))
                         + ( - qi(iV(i, j)) + qi(iV(i, j+1)) ) /dy;
         end
      
      %% Corners
      % bottom left
      i = 1;
      j = 1;
             qo(iP(i,j)) = (                + qi(iU(i+1, j)) ) /dx... % - qi(iU(i, j))
                         + (                + qi(iV(i, j+1)) ) /dy;  % - qi(iV(i, j))
      % bottom right
      i = nx;
      j = 1;
             qo(iP(i,j)) = ( - qi(iU(i, j))                  ) /dx... % + qi(iU(i+1, j))
                         + (                + qi(iV(i, j+1)) ) /dy;  % - qi(iV(i, j))

      % top left
      i = 1;
      j = ny;
             qo(iP(i,j)) = (                + qi(iU(i+1, j)) ) /dx... % - qi(iU(i, j))
                         + ( - qi(iV(i, j))                  ) /dy;  % + qi(iV(i, j+1))

      % top right
      i = nx;
      j = ny;
             qo(iP(i,j)) = ( - qi(iU(i, j))                  ) /dx... % + qi(iU(i+1, j))
                         + ( - qi(iV(i, j))                  ) /dy;  % + qi(iV(i, j+1))

end 


%% ##########	############################################################
function qo = OpGrad(qi)
    global np nu nx ny dx dy iU iV iP 
    %% Gradient operator
        % input:  P type (np elements) % 25 for 5*5 grid
        % output: u type (nu elments) % 40...

    %% Input size check:
    if size(qi, 2) ~= 1
        error('Grad:: input is not a column vector')
    end
    
    if size(qi, 1) ~= np
        error('Grad:: input size mismatch')
    end 
    
    %% Initialize output
    qo = nan(nu, 1);
    
    %% Inner domain
        %% X-direction gradient
    for i = 2:1:nx
        for j = 1:1:ny
            qo(iU(i,j)) = ( - qi(iP(i-1, j)) + qi(iP(i, j)) ) / dx ;
        end % say computing P-Grad at u22 location, we need P input at p12, p22
    end     % & there's no such thing as BC for P-input
    
        %% Y-direction gradient
    for i = 1:1:nx
        for j = 2:1:ny
            qo(iV(i,j)) = ( - qi(iP(i, j-1)) + qi(iP(i, j)) ) / dy ;
        end
    end
end


%% ##########	############################################################
function qo = OpLap(qi)
    global nu iU iV nx ny dx dy 
    %% Laplace operator
    % input: u-type (nu elements)
    % output: u-type (nu elements)
    
    %% Input size check:
    if size(qi, 2) ~= 1
        error('Lap:: input is not a column vector')
    end
    
    if size(qi, 1) ~=  nu
        error('Lap:: input size mismatch')
    end

    %% Initialize output
    qo = nan(nu, 1); % as output is u-type
    
  %% 1. u-component
    %% Inner Domain
    for i = 3:1:(nx-1)
        for j = 2:1:(ny-1) % y-wise, starting from 2 is ok
            qo(iU(i, j)) = ( qi(iU(i-1, j)) -2 *qi(iU(i, j)) + qi(iU(i+1, j)) ) /dx^2 ...
                         + ( qi(iU(i, j-1)) -2 *qi(iU(i, j)) + qi(iU(i, j+1)) ) /dy^2 ;
        end     % as dealing with x comp. there's only iU
    end

    %% Edges
    % left inner   (type 1)
    i = 2;  % remember, bcs 1st velocity nodes enters the 2nd cell
        for j = 2:1:(ny-1)
            qo(iU(i, j)) = (                 -2 *qi(iU(i, j)) + qi(iU(i+1, j)) ) /dx^2 ... % uBC_L
                         + ( +qi(iU(i, j-1)) -2 *qi(iU(i, j)) + qi(iU(i, j+1)) ) /dy^2 ;
        end
    
    % right inner   (type 1) 
    i = nx;
        for j = 2:1:(ny-1)
            qo(iU(i, j)) = ( +qi(iU(i-1, j)) -2 *qi(iU(i, j))                  ) /dx^2 ... % uBC_R
                         + ( +qi(iU(i, j-1)) -2 *qi(iU(i, j)) + qi(iU(i, j+1)) ) /dy^2 ;
        end
    
    % bottom inner  (type 2)
    j = 1; 
        for i = 3:1:(nx-1)      % remember the -ve symbol in 'bottom' position
            qo(iU(i, j)) = ( +qi(iU(i-1, j)) -2 *qi(iU(i, j)) + qi(iU(i+1, j)) ) /dx^2 ...
                         + ( -qi(iU(i, j))   -2 *qi(iU(i, j)) + qi(iU(i, j+1)) ) /dy^2 ; % 2*uBC_B
        end 
    
    % top inner     (type 2)
    j = ny;
        for i = 3:1:(nx-1)
            qo(iU(i, j)) = ( +qi(iU(i-1, j)) -2 *qi(iU(i, j)) + qi(iU(i+1, j)) ) /dx^2 ...
                         + ( +qi(iU(i, j-1)) -2 *qi(iU(i, j)) - qi(iU(i, j))   ) /dy^2 ; % 2*uBC_T
        end

    %% Corners
    % left bottom (t2 + t1)
    i = 2;
    j = 1;
            qo(iU(i, j)) = (                 -2 *qi(iU(i, j)) + qi(iU(i+1, j)) ) /dx^2 ... % uBC_L
                         + ( -qi(iU(i, j))   -2 *qi(iU(i, j)) + qi(iU(i, j+1)) ) /dy^2 ; % 2*uBC_B


    % right bottom (t2 + t1)
    i = nx;
    j = 1;
            qo(iU(i, j)) = ( +qi(iU(i-1, j)) -2 *qi(iU(i, j))                  ) /dx^2 ... % uBC_R
                         + ( -qi(iU(i, j  )) -2 *qi(iU(i, j)) + qi(iU(i, j+1)) ) /dy^2 ; % 2*uBC_B

    % left top (t2 + t1)
    i = 2;
    j = ny;
            qo(iU(i, j)) = (                 -2 *qi(iU(i, j)) + qi(iU(i+1, j)) ) /dx^2 ... % uBC_L
                         + ( +qi(iU(i, j-1)) -2 *qi(iU(i, j)) - qi(iU(i, j))   ) /dy^2 ; % 2*uBC_T

    % right top (t2 + t1)
    i = nx;
    j = ny;
            qo(iU(i, j)) = ( +qi(iU(i-1, j)) -2 *qi(iU(i, j))                  ) /dx^2 ... % uBC_R
                         + ( +qi(iU(i, j-1)) -2 *qi(iU(i, j)) - qi(iU(i, j))   ) /dy^2 ; % 2*uBC_T

  %% 2. v-Component
    %% Inner Domain
     for i = 2:1:(nx-1)
        for j = 3:1:(ny-1)
            qo(iV(i, j)) = ( qi(iV(i-1, j)) -2 *qi(iV(i, j)) + qi(iV(i+1, j)) ) /dx^2 ...
                         + ( qi(iV(i, j-1)) -2 *qi(iV(i, j)) + qi(iV(i, j+1)) ) /dy^2 ;
        end
     end

    %% Edges
    % left inner    (type 2)
    i = 1;
        for j = 3:1:(ny-1)    % watch
            qo(iV(i, j)) = ( -qi(iV(i, j))   -2 *qi(iV(i, j)) + qi(iV(i+1, j)) ) /dx^2 ...% 2*uBC_L
                         + ( +qi(iV(i, j-1)) -2 *qi(iV(i, j)) + qi(iV(i, j+1)) ) /dy^2 ;
        end
      
    % right inner  (type 2) 
    i = nx;
        for j = 3:1:(ny-1)                                      % watch
            qo(iV(i, j)) = ( +qi(iV(i-1, j)) -2 *qi(iV(i, j)) - qi(iV(i, j))   ) /dx^2 ...% 2*uBC_R
                         + ( +qi(iV(i, j-1)) -2 *qi(iV(i, j)) + qi(iV(i, j+1)) ) /dy^2 ;
        end
 
    % bottom inner  (type 1)
    j = 2;
        for i = 2:1:(nx-1)
            qo(iV(i, j)) = ( +qi(iV(i-1, j)) -2 *qi(iV(i, j)) + qi(iV(i+1, j)) ) /dx^2 ...
                         + (                 -2 *qi(iV(i, j)) + qi(iV(i, j+1)) ) /dy^2 ; % uBC_B
        end
 
    % top inner     (type 1)
    j = ny;
        for i = 2:1:(nx-1)
            qo(iV(i, j)) = ( +qi(iV(i-1, j)) -2 *qi(iV(i, j)) + qi(iV(i+1, j)) ) /dx^2 ...
                         + ( +qi(iV(i, j-1)) -2 *qi(iV(i, j))                  ) /dy^2 ; % uBC_T
        end

    %% Corners
    % left bottom
    i = 1;
    j = 2;                   % watch
            qo(iV(i, j)) = ( -qi(iV(i, j))   -2 *qi(iV(i, j)) + qi(iV(i+1, j)) ) /dx^2 ...% 2*uBC_L
                         + (                 -2 *qi(iV(i, j)) + qi(iV(i, j+1)) ) /dy^2 ; % uBC_B
    
    % right bottom
    i = nx;
    j = 2;                                                    % watch
            qo(iV(i, j)) = ( +qi(iV(i-1, j)) -2 *qi(iV(i, j)) - qi(iV(i, j))   ) /dx^2 ...% 2*uBC_R
                         + (                 -2 *qi(iV(i, j)) + qi(iV(i, j+1)) ) /dy^2 ; % uBC_B

    % left top
    i = 1;
    j = ny;                   % watch
            qo(iV(i, j)) = ( -qi(iV(i, j))   -2 *qi(iV(i, j)) + qi(iV(i+1, j))  ) /dx^2 ...% 2*uBC_L
                         + ( +qi(iV(i, j-1)) -2 *qi(iV(i, j))                   ) /dy^2 ; % uBC_T
    
    % right top
    i = nx;
    j = ny;                                                   % watch
            qo(iV(i, j)) = ( +qi(iV(i-1, j)) -2 *qi(iV(i, j)) - qi(iV(i, j))   ) /dx^2 ...% 2*uBC_R
                         + ( +qi(iV(i, j-1)) -2 *qi(iV(i, j))                  ) /dy^2 ; % uBC_T

end

%% ##########	############################################################
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


%% ##########	############################################################
function [Y, U] = PlotGenNmr(Vec, fn)
    % first, clips the solution into two (u, v), then form matrix
    % clips a vector in the midpoint of matrix

    global nx ny Lx Ly...
    uBC_L uBC_R uBC_B uBC_T vBC_L vBC_R vBC_T vBC_B
    
    [Mat_u, Mat_v] = VectToMat(Vec);
        % clipping the input vect' into 2, n turning them into mat'
    
    Mat_u(1,:) = uBC_L; % Replacing the BCs within sol'
%     Mat_v(:,1) = vBC_B;
        % This si redundant part thought

    Mat_u_vw = flipud(Mat_u'); % making similar to diagram
%     Mat_v_vw = flipud(Mat_v');
    
    %% Clipping-off u1 values
        % ..........across channel (along x2 line), at channel center
   if mod(nx,2) == 1   % odd num'
       q = (nx+1)/2;    % middle, not CG line
       K1 = Mat_u_vw(:,q);
       K2 = Mat_u_vw(:,q+1); % clipping off the col' across CG
       K_vw = (K1 + K2)/2;
    
    else                % even num'
       q = nx /2;       % middle, not CG line
       K_vw = Mat_u_vw(:, q+1); % clipping off col' across CG
    end
    
    %% Plotter
    u1_y_nm = flipud(K_vw); % reshaping back, for plotting
    u1_y_nm_ex = [uBC_B(q); u1_y_nm; uBC_T(q)];
%     u1_y_nm(ny+1) = uBC_T(q);   % appending the boundary value
    
        %% Grid
    y0 = linspace(-Ly/2, Ly/2, 2*ny +1); % additional 1 for boundaries
%     y_ov = y0(1:2:end-1); % excluding the top boundary
    y_eu = y0(2:2:end); % takes up even values
    y_eu_ex = [-Ly/2, y_eu, Ly/2]; % appended boundary nodes
        
        %% The Plot
    figure(fn)
    plot(y_eu_ex, u1_y_nm_ex, 'LineWidth', 1.5), grid on
    xlabel('Channel Width, $$ x_2(m) $$', 'Interpreter', 'latex')
    ylabel('Velocity, $$ u_1(m/s) $$', 'Interpreter', 'latex')
    title('$$ x_2 $$ vs $$ u_1 $$', 'Interpreter', 'latex')
    
    %% Output
    Y = y_eu_ex; U = u1_y_nm_ex;

end

%% ##########	############################################################
function [X] = PrintMat(Mat, fn,  chtit, chrx, chry)
    % Mat = input Matrix; nx, ny = matrix size; fn = figure no.; chr input

    Print = Mat';  % transpose the matrix to make it similar to diagram
    
    Mat_vw = flipud(Mat'); % Making'em like diagram
    
    X = Mat_vw;

    %% Prints the u(x) components
    figure(fn)
    contourf(Print, 14, 'ShowText','off' )
    title([num2str(chtit)], 'Interpreter', 'latex')
    xlabel([num2str(chrx)], 'Interpreter', 'latex')
    ylabel([num2str(chry)], 'Interpreter', 'latex')

end


%% ##########	############################################################
function [Vec_u, Vec_v] = VectClipper(Vec, nx, ny)
    % clips the input vector into u and v component
    nu = (nx-1) *ny + nx *(ny-1);

 %% Clipping the long vector 
    str_u = 1;              end_u = (nx-1)*ny;
    str_v = (nx-1)*ny +1;   end_v = nu;
    Vec_u = Vec(str_u :end_u);
    Vec_v = Vec(str_v :end_v);
end


%% ##########	############################################################
function [Mat_u, Mat_v] = VectToMat(Vec)
    % transform, vector to matrix
    % keep it into its original form, not transpose to look-like diag'
global nx ny nu
    nu = (nx-1) *ny + nx *(ny-1);

 %% Clipping the long vector 
    str_u = 1;              end_u = (nx-1)*ny;
    str_v = (nx-1)*ny +1;   end_v = nu;
    Vec_u0 = Vec(str_u :end_u);
    Vec_v0 = Vec(str_v :end_v);

  %% Memory Allocation
    iU = nan(nx, ny);
    iV = nan(nx, ny);
    
    %% Long vector to nx*ny Matrix
        % with boundary values shown
    id_u = 1;
    for i = 2:1:nx
        for j = 1:1:ny
            iU(i,j) = Vec_u0(id_u);
            id_u = id_u +1;
        end 
    end
    Mat_u = iU; % output 
    
    id_v = 1;
    for i = 1:1:nx
        for j = 2:1:ny
            iV(i,j) = Vec_v0(id_v);
            id_v = id_v +1;
        end 
    end
    Mat_v = iV; % output
end


%% ##########	############################################################
function [Mat_p] = VectToMatP(Vec)
    % transform, vector to matrix
    % keep it into its original form
global nx ny
    
  %% Memory Allocation
    iP = nan(nx, ny);
    
    %% Long vector to nx*ny Matrix
        % with boundary values shown
    id_p = 1;
    for i = 1:1:nx
        for j = 1:1:ny
            iP(i,j) = Vec(id_p);
            id_p = id_p +1;
        end 
    end
    Mat_p = iP; % output

end
