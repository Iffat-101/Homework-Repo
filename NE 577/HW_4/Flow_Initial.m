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
            Uu_o1_Mt(i,j) = 0.0125; % no pre-u in the domain
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