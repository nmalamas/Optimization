function [f_out] = f(x, P, q)
% The quadratic function we study

    f_out = (1/2)*x.'*P*x + q.'*x;

end