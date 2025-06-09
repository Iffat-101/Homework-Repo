function y = func_intrp(x1, x2, y1, y2, x)
% this is a 2 point interpolation function
% takes up values of 2 points at returns value for x
    y = (x-x1)/(x1-x2) *(y1-y2) + y1;
end 


