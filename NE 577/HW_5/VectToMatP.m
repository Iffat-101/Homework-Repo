function [Mat_p] = VectToMatP(Vec)
    % transform, vector to matrix
    % keep it into its original form
global nx ny
    
  %% Memory Allocation
    iP = nan(nx, ny);
    
    %% Long vector to nx*ny Matrix
        % with boundary values shown
    id_p = 1;
    for i = 1:1:nx
        for j = 1:1:ny
            iP(i,j) = Vec(id_p);
            id_p = id_p +1;
        end 
    end
    Mat_p = iP; % output

end