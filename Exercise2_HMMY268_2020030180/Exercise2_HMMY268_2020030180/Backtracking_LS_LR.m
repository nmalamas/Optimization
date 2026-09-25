function [theta_opt] = Backtracking_LS_LR(X_augm, y, theta_0, lambda, epsilon, alpha, beta, n, theta_star)

theta = theta_0;
grad = grad_J_R(X_augm, theta, lambda, y, n);
fun_val = J_R(X_augm, theta, lambda, y, n);

iters = 0;

while (norm(grad) > epsilon)
    iters = iters + 1;
    t = 1;

    while (fun_val - J_R(X_augm, theta - t*grad, lambda, y, n) < alpha*t*norm(grad)^2)
        t = beta*t;
    end

    theta = theta - t*grad;
    fun_val = J_R(X_augm, theta, lambda, y, n);
    grad = grad_J_R(X_augm, theta, lambda, y, n);

    fprintf('Relative Percentage Error: %.6f [%%]', ( norm(theta - theta_star)/norm(theta_star) ) );
    fprintf('\n');
end

theta_opt = theta;
if (theta_star - epsilon <= theta_opt & theta_opt <= theta_star + epsilon)
    fprintf('\n\n');
    fprintf('The GD algorithm with BT-LS found the solution of the LR problem to be the same as fmincon()!\n');
    fprintf('θ_* = [%.3f %.3f %.3f]\n', theta_opt(1), theta_opt(2), theta_opt(3));
end

end