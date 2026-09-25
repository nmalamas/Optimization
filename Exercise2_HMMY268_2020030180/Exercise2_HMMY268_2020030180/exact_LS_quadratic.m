function [xk,f_xk] = exact_LS_quadratic(P, q, x0, epsilon)
    
    x = x0;
    xk = [x];
    f_xk = [f(x0, P, q)];
    iters = 0;
    
    grad = g(x0, P, q);
    
    while (norm(grad) > epsilon)
        
        % Find t_* using the closed form
        t = (norm(grad)^2)/(grad.'*P*grad); 

        % Update
        x = x - t*grad;
        xk = [xk x];
        grad = g(x, P, q);       
        f_xk = [f_xk f(x, P, q)];
        
        iters = iters + 1; % increase iteration counter
    end

end