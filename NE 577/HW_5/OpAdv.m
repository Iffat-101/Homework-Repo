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
