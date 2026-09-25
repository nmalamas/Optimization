function [h_val] = h(A,b,x)

    m = length(A(:,1));
    sum = 0;

    for i = 1:m
        % Get the rows of A
        a_i_T = A(i,:);
        scal = 1/(b(i) - a_i_T * x)^2;

        sum = sum + ( scal * a_i_T.' * a_i_T );
    end

    h_val = sum;

end