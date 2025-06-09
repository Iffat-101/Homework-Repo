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
