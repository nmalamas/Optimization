function grad = gradient_SVM_barrier( w_inner_iter, X_augm, y, t )

N = length(y);
sum = 0;

for i = 1:N
    scal = y(i) / (1 - y(i) .* X_augm(:,i)' * w_inner_iter);
    sum = sum + scal .* X_augm(:,i);
end

grad = t*w_inner_iter + sum;

end