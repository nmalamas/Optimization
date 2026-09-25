% Convex Optimization 2025-26 %
% --- EXERCISE 1 --- %
% Nikolaos Malamas (2020030180)

%%
%--- FIRST PART ---%
% A. Done in paper - results in the report

%%
clc; close all; clear all;
fprintf('\t\t### First Part ###\n');

% B.
% (a)
n = 2;

% i.(a)
A = randn(n,n);

[U,S,V] = svd(A);

UU_trans = U*U.';
for i=1:n
    for j=1:n
        UU_trans(i,j) = round(UU_trans(i,j),3);
    end
end

fprintf('+==========================================+\n');
if (UU_trans == eye(n))
    fprintf("U*U' is equal to the identity matrix!\n");
end

fprintf('+==========================================+\n');
U_transU = U.'*U;
for i=1:n
    for j=1:n
        U_transU(i,j) = round(U_transU(i,j),3);
    end
end

if (U_transU == eye(n))
    fprintf("U'*U is equal to the identity matrix!\n");
end
fprintf('+==========================================+\n');

% i.(b)
lambda_min = 1;
lambda_max = 100;

% Generate a range of eigenvalues for the optimization problem
z = lambda_min + (lambda_max - lambda_min)*rand(n-2,1);

% Generate the vector of the n eigenvalues
eig_P = [lambda_min; lambda_max; z];

% Construct matrix Lambda
Lambda = diag(eig_P);

% Define the condition number of the problem
K = lambda_max/lambda_min;

% ii.
% Construct random vector q
q = rand(n,1);

% Construct random positive definite matrix P
% using condition number K
P = U*Lambda*U.';

% iii.
% Find x_* using the closed form solution
x_star = - inv(P) * q;

% Calculate the optimal value, based on the closed for solution
f_star = f(x_star, P, q); 

% iv.
A = [];
b = [];
x0 = rand(n,1);

% Use fmincon to minimize the function f,
% with no constraints
x_mincon = fmincon(@(x)(1/2)*x.'*P*x + q.'*x, x0, A, b);

fprintf('+==========================================+\n');
if(round(x_star,4) == round(x_mincon,4))
   fprintf('fmincon gives the same solution as the closed form solution!\n');
end
fprintf('+==========================================+\n');

% v.
epsilon = 10^(-6);
x0 = 5*ones(n,1); % Woud work even for x0 = randn(n,1),
% but with less noticable visual effects when plotting the level sets

% GD with exact line search
[X_exact, F_exact] = exact_LS_quadratic(P, q, x0, epsilon);

% GD w/Backtracking line search
alpha = 0.1;
beta = 0.7;

[X_bt, F_bt] = Backtracking_LS_quadratic(P, q, x0, epsilon, alpha, beta);

% vi.
if (n == 2)
    step = 0.1;
    x1 = min(x0(1), x_star(1)) - 5:step:max(x0(1), x_star(1)) + 5;
    x2 = min(x0(2), x_star(2)) - 5:step:max(x0(2), x_star(2)) + 5;
    for i1 = 1:length(x1)
        for i2 = 1:length(x2)
            x_tmp = [x1(i1); x2(i2)];
            f_(i1,i2) = f(x_tmp, P, q);
        end
    end 
    
    figure
    if (length(F_exact) < 20)
        ls = linspace(F_exact(1),F_exact(end),10);
        contour(x1, x2, f_', ls);
        axis('square'); % plot the level sets of f_duadratic at the positions x_k
        hold on   
        plot(X_exact(1,:), X_exact(2,:), '-or'); % plot the trajectory of x_k
        hold off
    else
        contour(x1, x2, f_', sort(F_exact));
        axis('square'); % plot the level sets of f_duadratic at the positions x_k
        hold on   
        plot(X_exact(1,:), X_exact(2,:), '-or'); % plot the trajectory of x_k
        hold off
    end
      
    figure
    if (length(F_bt) == 2)  
        ls = linspace(F_exact(1),F_exact(end),10);
        contour(x1, x2, f_', ls);
        axis('square'); % plot the level sets of f_duadratic at the positions x_k
        hold on   
        plot(X_bt(1,:), X_bt(2,:), '-or'); % plot the trajectory of x_k
        hold off
    else
        contour(x1, x2, f_', sort(F_bt));
        axis('square'); % plot the level sets of f_duadratic at the positions x_k
        hold on   
        plot(X_bt(1,:), X_bt(2,:), '-or'); % plot the trajectory of x_k
        hold off
    end
end

% vii.
figure
semilogy([1:length(F_exact)]', (F_exact-F_exact(end)), 'b', 'LineWidth', 1.5);
hold on;
semilogy([1:length(F_bt)]', (F_bt-F_bt(end)), 'r', 'LineWidth', 1.5);
grid on;
xlabel('$k$', 'FontSize', 13, 'Interpreter','latex');
ylabel('$log(f(\mathbf{x}_k)- p_*)$', 'FontSize', 13, 'Interpreter', 'latex');
legend('Exact LS', 'Backtracking LS');
axis tight

% Calculate the slope of the plot
c = 1 - (1/K);
slope = log10(c);
fprintf('The slope of the plot for K=%d is: %.3f\n', K, slope);
fprintf('+==========================================+\n');

% viii.
k_epsilon = K*log((f(x0,P,q) - f_star) / epsilon);
fprintf('Max #iters needed for 10^(%d)-close convergence, with K=%d, is: %.2f\n', log10(epsilon), K, k_epsilon);
fprintf('+==========================================+\n');
fprintf('\n');

%%
%--- SECOND PART ---%
fprintf('\t\t### Second Part ###');

% A. Data Generation
N = 2; % #dimensions_of_datapoints
n = 100; % #datapoints

% 1. 
w = rand(N, 1);
b = rand;

% Un-comment to get varying values of std
% for std = 0.1:0.4:0.9 
for std = 0.15:0.15
    
    % 2.
    % x_0 exists in the Hyperplane: H = {w^T * x = b}
    x0 = (w/norm(w).^2) .* b;
    
    % 3.
    % Randomly generate the {0,1} class labels
    y = zeros(1,n);
    for i=1:n
        y(i) = round(rand);
    end
    
    % 4.
    % Randomly generate the data points in the N-dimensional space
    X = zeros(N,n);
    
    for i=1:n
        if (y(i) == 1)
            X(:,i) = x0 + w + std * randn(N,1);
        else
            X(:,i) = x0 - w + std * randn(N,1);
        end
    end
    
    % 5.
    step = 0.01;
    x1 = x0(1) - 5: step: x0(1) + 5-step;
    x2 = x0(2) - 5: step: x0(2) + 5-step;
    x = [x1; x2];
    
    [a,c] = line_func_2D(w, b, x0);
    sep_line = a*x + c;
    
    if (N == 2)        
        % Plot the data points and the separating hyperplane
        figure;
        scatter(X(1, y == 0), X(2, y == 0), 'r', 'filled'); % Class 0
        hold on;
        scatter(X(1, y == 1), X(2, y == 1), 'b', 'filled'); % Class 1
        hold on;
    
        plot(x(1:10:end), sep_line(1:10:end), 'Color', 'k', 'LineStyle', ':');
        hold on;
        plot(x0(1), x0(2), 'Marker', 'o', 'MarkerEdgeColor', 'm', 'MarkerFaceColor', 'm', 'MarkerSize', 5);
        text(x0(1) + 0.01, x0(2) + 0.01, 'x_0', 'FontSize', 12 ,'Color', 'k', 'FontWeight', 'bold');
    
        xlabel('x_1','FontSize',13);
        ylabel('x_2','FontSize',13);
        x1_min = min(X(1,:))-0.2; x1_max = max(X(1,:))+0.2; x2_min = min(X(2,:))-0.2; x2_max = max(X(2,:))+0.2;
        axis([x1_min x1_max x2_min x2_max]);
        lgd = legend('Class 0', 'Class 1');
        title(lgd, 'True Labels');
        grid on;
        hold off;
        title(['Seperable Dataset with $\sigma^2=$', num2str(std)] ,'Interpreter','latex');
    end
    
    % 6.
    % Define vector θ := (b,w)
    theta_0 = [b; w];
    
    % Augment dataset X 
    X_augm = zeros(N+1,n);
    
    for i=1:n
        X_augm(:,i) = [-1; X(:,i)];
    end
end

%%
% B. Classification via Logistic Regression

% 1.
% Iterate through (positive) values of parameter lambda
lambda = 10^0;

% Un-comment to get results for varying values of lambda
%while (lambda > 10^(-2))
    
    
    % Define null constraints for fmincon
    A = [];
    b = [];
    
    % Define the initial theta of the LR
    theta_star = fmincon(@(theta) (-1/n)*sum( y*log(1./(1+exp(-X_augm.'*theta))) ...
                + (1-y)*log(1-(1./(1+exp(-X_augm.'*theta)))) ) ...
                + lambda .* (1/2)*norm(theta)^2, theta_0, A, b);
    
    
    if (N == 2)
        b_star = theta_star(1);
        w_star = theta_star(2:3);
        
        [a,c] = line_func_2D(w_star, b_star, x0);
        sep_line_opt = a*x + c;
        
        % Plot the data points and the separating hyperplane
        figure;
        scatter(X(1, y == 0), X(2, y == 0), 'r', 'filled'); % Class 0
        hold on;
        scatter(X(1, y == 1), X(2, y == 1), 'b', 'filled'); % Class 1
        hold on;
    
        plot(x(1:10:end), sep_line(1:10:end), 'Color', 'k', 'LineStyle', ':');
        hold on;
        plot(x(1:10:end), sep_line_opt(1:10:end), 'Color', 'k', 'LineStyle', ':', 'LineWidth', 1.7);
        hold on;
        plot(x0(1), x0(2), 'Marker', 'o', 'MarkerEdgeColor', 'm', 'MarkerFaceColor', 'm', 'MarkerSize', 5);
        text(x0(1) + 0.01, x0(2) + 0.01, 'x_0', 'FontSize', 12 ,'Color', 'k', 'FontWeight', 'bold');
        lgd1 = legend('Initial seperating line', 'Optimal seperating line');
        hold on
    
        xlabel('x_1','FontSize',13);
        ylabel('x_2','FontSize',13);
        x1_min = min(X(1,:))-0.2; x1_max = max(X(1,:))+0.2; x2_min = min(X(2,:))-0.2; x2_max = max(X(2,:))+0.2;
        axis([x1_min x1_max x2_min x2_max]);
        lgd2 = legend('Class 0', 'Class 1');
        title(lgd2, 'True Labels');
        grid on;
        hold off;
        title(['Seperable Dataset with $\sigma^2=$', num2str(std), ' and $\lambda=$', num2str(lambda)] ,'Interpreter','latex');
    end

 %lambda = 0.1*lambda;
 %end


% 2.
fprintf('+==========================================+\n');
fprintf('BackTracking LS applied for the optimal of J_R(θ)\n');
fprintf('-------------------------------------------\n');
theta_opt = Backtracking_LS_LR(X_augm, y, theta_0, lambda, epsilon, alpha, beta, n, theta_star);
fprintf('+==========================================+\n');
