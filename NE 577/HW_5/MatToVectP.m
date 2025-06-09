function Vec_p = MatToVectP(Mat_p)
    global nx ny nu np 
        % here, boundary values are 'Not' clipped off for u and v
        % then Matrix is converted into row vector
    %% Memory Allocation
    V_pL = nan(1, np);

    %% Creating Long vector
    id_p = 1;
    for i = 1:1:nx  % retains 1st node values
        for j = 1:1:ny
            V_pL(id_p) = Mat_p(i,j);
            id_p = id_p +1;
        end
    end
    Vec_p = V_pL;

end 
 