function [iP, iU, iV] = GenPointer(nx, ny)
% global nx ny
%% Memory Allocation
    iP =  nan(nx, ny);
    iU =  nan(nx, ny); 
    iV =  nan(nx, ny);  % allocating NaN for variables
    
    %% Pointer matrix for P
    id_p = 1; % index to be used in vector variable P
    for i =1:1:nx
        for j = 1:1:ny
            iP(i,j) = id_p; % putting number in pointers
            id_p = id_p + 1;       
        end
    end 
%     iP; % output

    %% Pointer Matrix for ux
    id_u = 1; % index to be used in vector variable u = [ux; uy]
    for i =2:1:nx   % starts from 2
        for j = 1:1:ny
            iU(i,j) = id_u;
            id_u = id_u + 1;
        end
    end
%     iU; % output

    %% Pointer Matrix for vy
    for i =1:1:nx
        for j = 2:1:ny  % starts from 2
            iV(i,j) = id_u;
            id_u = id_u + 1;   % long vector index    
        end
    end
%     iV; % output

    %% Visualizer
    iP_vw = flipud(iP'); % viewing the memory locations as like diagram
    iU_vw = flipud(iU');
    iV_vw = flipud(iV');
end
