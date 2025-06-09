function [X, Y] = BCAppend(Vec)
    % Vec = input Vector; nx, ny = matrix size; fn = figure no.; chr input
    % first, clip the vector into two; generate matrix > contour for 'em
    
    global uBC_L uBC_R uBC_B uBC_T vBC_L vBC_R vBC_T vBC_B
    [Mat_u, Mat_v] = VectToMat(Vec);
    
    Mat_u(1,:) = uBC_L; % Replacing the BCs
    Mat_v(:,1) = vBC_B;

    Mat_uT = Mat_u';  % transpose the Matrix makin it similar to diagrm
    Mat_uT = [Mat_uT, uBC_R']; % appending right BC
    uBC_Te = [uBC_T, 0]; % adding extra 0 to matchup size
    uBC_Be = [uBC_B, 0]; 
    Mat_uT = [uBC_Be; Mat_uT; uBC_Te]; % Appending the boundaries
        % these are half-grid appendix, so, makes domain look longer
    
    Mat_vT = Mat_v';
    Mat_vT = [Mat_vT; vBC_T];
    vBC_Le = [vBC_L, 0];
    vBC_Re = [vBC_R, 0];
    Mat_vT = [vBC_Le', Mat_vT, vBC_Re'];

    X = Mat_uT'; % reverting back to Matlab form
    Y = Mat_vT';

    Vw_1 = flipud(Mat_uT); % for viewing how it looks like in mat form
    Vw_2 = flipud(Mat_vT);
    
end
