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
