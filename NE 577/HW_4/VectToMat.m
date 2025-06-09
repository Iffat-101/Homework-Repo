function [Mat_u, Mat_v] = VectToMat(Vec)
    % transform, vector to matrix
    % keep it into its original form, not transpose to look-like diag'
global nx ny nu
    nu = (nx-1) *ny + nx *(ny-1);

 %% Clipping the long vector 
    str_u = 1;              end_u = (nx-1)*ny;
    str_v = (nx-1)*ny +1;   end_v = nu;
    Vec_u0 = Vec(str_u :end_u);
    Vec_v0 = Vec(str_v :end_v);

  %% Memory Allocation
    iU = nan(nx, ny);
    iV = nan(nx, ny);
    
    %% Long vector to nx*ny Matrix
        % with boundary values shown
    id_u = 1;
    for i = 2:1:nx
        for j = 1:1:ny
            iU(i,j) = Vec_u0(id_u);
            id_u = id_u +1;
        end 
    end
    Mat_u = iU; % output 
    
    id_v = 1;
    for i = 1:1:nx
        for j = 2:1:ny
            iV(i,j) = Vec_v0(id_v);
            id_v = id_v +1;
        end 
    end
    Mat_v = iV; % output
end