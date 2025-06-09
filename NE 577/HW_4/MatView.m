function [Mat_u_vw, Mat_v_vw] = MatView(Mat_u, Mat_v)
    % transform the vectors into how they look like in the diagram
    Mat_u_vw = flipud(Mat_u');
    Mat_v_vw = flipud(Mat_v');
end