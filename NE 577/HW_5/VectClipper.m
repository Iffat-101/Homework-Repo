function [Vec_u, Vec_v] = VectClipper(Vec, nx, ny)
    % clips the input vector into u and v component
    nu = (nx-1) *ny + nx *(ny-1);

 %% Clipping the long vector 
    str_u = 1;              end_u = (nx-1)*ny;
    str_v = (nx-1)*ny +1;   end_v = nu;
    Vec_u = Vec(str_u :end_u);
    Vec_v = Vec(str_v :end_v);
end

