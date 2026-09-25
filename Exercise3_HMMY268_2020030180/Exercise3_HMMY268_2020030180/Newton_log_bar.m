function [xk,f_xk] = Newton_log_bar(A, b, c, x0, epsilon, alpha, beta)

    x = x0;
    xk = [x];

    [~,~,fun_val] = f(A, b, c, x0);
    f_xk = [fun_val];
    grad = g(A, b, c, x0);
    hess = h(A,b,x0);

    iters = 0;
 
    while (1)
        iters = iters + 1;
        fprintf('Iteration: %d\n', iters);
        
        Dx_Nt = -inv(hess) * grad;
        lambda_sq = grad.' * inv(hess) * grad;

        if (lambda_sq/2 <= epsilon)
            break;
        end
        
        t = 1;
        % Iteratively search for a point in the domain of f
        while ( any(b - A*(x + t*Dx_Nt) <= 0) )
            x_tmp = x + t*Dx_Nt;
            fprintf('x_k+1 = [%.4f,%.4f] does NOT belong to dom{f} \n', x_tmp(1), x_tmp(2));
    
            t = beta * t;
        end

        fprintf('----------------\n');
    
        % After finding a valid point, proceed with the standard BT-LS
        [~,~,f0] = f(A, b, c, x + t*Dx_Nt);
        while (f0 > fun_val +  alpha*t*norm(grad)^2)
            t = beta * t;
            [~,~,f0] = f(A, b, c, x + t*Dx_Nt);
        end

        % Update
        x = x + t*Dx_Nt;
        xk = [xk x];
        [~,~,fun_val] = f(A, b, c, x);
        f_xk = [f_xk fun_val];
        grad = g(A, b, c, x);
        hess = h(A,b,x);

    end

end