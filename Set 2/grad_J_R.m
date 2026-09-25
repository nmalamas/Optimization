function [grad_J_R_out] = grad_J_R(X_augm, theta, lambda, y, n)

% Gradient of J(theta)
grad_J = [];
% for i=1:length(X_augm(1,:))
%     grad_J = grad_J + (h_theta(X_augm(:,i), theta) - y(i))*X_augm(:,i);
% end
% grad_J = (1/n)*grad_J;
grad_J = (1/n)*X_augm*( 1./(1 + exp(-X_augm.'*theta)) - y.' );

% Gradient of the squared Euclidean norm R(theta)
grad_R = lambda * theta;

% Total gradient of l2-reguralized LR
grad_J_R_out = grad_J + grad_R;

end