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
