function [xk,f_xk] = Backtracking_LS_quadratic(P, q, x0, epsilon, alpha, beta)

    x = x0;
    xk = [x];

    fun_val = f(x0, P, q);
    f_xk = [fun_val];
    grad = g(x0, P, q);

    iters = 0;
 
    while (norm(grad) > epsilon)
        iters = iters + 1;
        t = 1;

        while (fun_val - f(x - t*grad, P, q)  <  alpha*t*norm(grad)^2)
            t = beta * t;
        end
    
        % Update
        x = x - t*grad;
        xk = [xk x];

        fun_val = f(x, P, q);
        f_xk = [f_xk fun_val];
        grad = g(x, P, q);

    end

end