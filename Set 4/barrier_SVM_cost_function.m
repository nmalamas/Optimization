function g_val = barrier_SVM_cost_function( w_new, X_augm, y, t )

f_0 = (1/2).*norm(w_new)^2;

phi = -sum(log(y' .* (X_augm' * w_new) - 1));

g_val = t*f_0 + phi;

end