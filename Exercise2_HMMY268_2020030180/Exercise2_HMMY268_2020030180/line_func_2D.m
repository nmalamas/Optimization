function [slope,c] = line_func_2D(w, b, x0)
% Trying to get to the form y = a*x + c

% Find the slope by getting 2 points on each axis: (0,b/w2) and (b/w1,0)
p1 = [0; b/w(2)];
p2 = [b/w(1); 0];

slope = (p1(2) - p2(2)) / (p1(1) - p2(1));

% Find c 
c = x0(2) - slope * x0(1);

end