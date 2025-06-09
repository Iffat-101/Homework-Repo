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
