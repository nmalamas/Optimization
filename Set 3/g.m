function [g_val] = g(A, b, c, x)

    m = length(A(:,1));
    sum = 0;

    for i = 1:m
        % Get the rows of A
        a_i_T = A(i,:);
        scal = 1/(b(i) - a_i_T * x);

        sum = sum + ( scal * a_i_T.' );
    end

    g_val = c + sum;

end