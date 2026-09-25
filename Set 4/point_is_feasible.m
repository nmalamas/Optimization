function flag = point_is_feasible(w_init, X_augm, y)
flag = 1;

N = length(y);

tmp = y' .* (X_augm' * w_init) > 1;
if ( sum(tmp) < N )
    flag = 0;
end

end