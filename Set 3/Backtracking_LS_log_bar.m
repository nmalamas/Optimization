function [xk,f_xk] = Backtracking_LS_log_bar(A, b, c, x0, epsilon, alpha, beta, x_opt)

    x = x0;
    xk = [x];

    [~,~,fun_val] = f(A, b, c, x0);
    f_xk = [fun_val];
    grad = g(A, b, c, x0);

    iters = 0;
 
    while (norm(x - x_opt) > epsilon)
        iters = iters + 1;
        fprintf('Iteration: %d\n', iters);
        t = 1;


        % Iteratively search for a point in the domain of f
        while ( any(b - A*(x - t*grad) <= 0) )
            x_tmp = x - t*grad;
            fprintf('x_k+1 = [%.4f,%.4f] does NOT belong to dom{f} \n', x_tmp(1), x_tmp(2));
    
            t = beta * t;
        end
    
        % After finding a valid point, proceed with the standard BT-LS
        [~,~,f0] = f(A, b, c, x - t*grad);
        while (fun_val - f0 <  alpha*t*norm(grad)^2)
            t = beta * t;
            [~,~,f0] = f(A, b, c, x - t*grad);
        end

        fprintf('----------------\n');

        % Update
        x = x - t*grad;
        xk = [xk x];
        [~,~,fun_val] = f(A, b, c, x);
        f_xk = [f_xk fun_val];
        grad = g(A, b, c, x);

    end

end