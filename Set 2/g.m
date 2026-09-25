function [g_out] = g(x, P, q)
% The gradient of the quadratic function we study

    g_out = P*x + q;

end