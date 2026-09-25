function [x1_val,x2_val,f_val] = f(A, b, c, x)

    infty = 10^3;
    NID = -1; % Not In Domain 

    if (b - A*x > zeros(length(b),1))
        f_val = c.' * x - sum(log(b - A*x));
        x1_val = x(1);
        x2_val = x(2);
    else
        f_val = infty;
        x1_val = NID;
        x2_val = NID;
    end

end