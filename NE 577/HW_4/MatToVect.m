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
