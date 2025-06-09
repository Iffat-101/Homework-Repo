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