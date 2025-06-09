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