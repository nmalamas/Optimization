function [J_R_out] = J_R(X_augm, theta, lambda, y, n)

J_R_out =  (-1/n)*sum( y*log(1./(1+exp(-X_augm().'*theta)) ) ...
           + (1-y)*log(1-(1./(1+exp(-X_augm.'*theta)))) ) ...
           + lambda .* (1/2)*norm(theta)^2;

end