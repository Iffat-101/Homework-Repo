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