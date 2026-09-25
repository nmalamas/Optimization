function hess = Hess_SVM_barrier( w_inner_iter, X_augm, y, t )

N = length(y);
sum = 0;

for i = 1:N
    scal = 1 / (1 - y(i) .* X_augm(:,i)' * w_inner_iter)^2;
    sum = sum + scal .* X_augm(:,i) * X_augm(:,i)';
end

hess = t*ones(length(X_augm(:,1))) + sum;

end