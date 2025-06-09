function [X] = PrintMat(Mat, fn,  chtit, chrx, chry)
    % Mat = input Matrix; nx, ny = matrix size; fn = figure no.; chr input

    Print = Mat';  % transpose the matrix to make it similar to diagram
    
    Mat_vw = flipud(Mat'); % Making'em like diagram
    
    X = Mat_vw;

    %% Prints the u(x) components
    figure(fn)
    contourf(Print, 14, 'ShowText','off' )
    title([num2str(chtit)], 'Interpreter', 'latex')
    xlabel([num2str(chrx)], 'Interpreter', 'latex')
    ylabel([num2str(chry)], 'Interpreter', 'latex')

end
