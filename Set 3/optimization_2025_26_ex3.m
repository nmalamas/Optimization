% Convex Optimization 2025-26 %
% --- EXERCISE 3 --- %
% Nikolaos Malamas (2020030180)

clc; clear all; close all;

% FIRST PART %
n = 2;
m = 20;

% 1. Create A, b, and c randomly
A = randn(m,n);
b = rand(m,1);
c = randn(n,1);

% 2. Minimize f using cvx
cvx_begin
    variables x_cvx(n)
    minimize(c' * x_cvx - sum(log(b - A * x_cvx)))
cvx_end

fprintf('Optimal point as given by CVX (x_cvx): (%.4f, %.4f)\n', x_cvx(1), x_cvx(2));
fprintf('\n');

% 3. Plot f and its level sets near the optimum point
if (n == 2)    
    step = 10^(-3); % 10^-4 takes up to 3 minutes to run
    x1 = x_cvx(1) - 0.25: step: x_cvx(1) + 0.25 - step;
    x2 = x_cvx(2) - 0.25: step: x_cvx(2) + 0.25 - step;
    [X1,X2] = meshgrid(x1, x2);

    % Calculate ALL the values of f,
    % even in points outside the domain
    for i = 1:length(X1(:,1))
        for j = 1:length(X1(:,1))
            x_tmp = [X1(i,j); X2(i,j)];
            [X1_(i,j), X2_(i,j), f_(i,j)] = f(A, b, c, x_tmp);
        end
    end   

    % Full version of f
    f_log_bar = f_;

    % Exclude the points that do NOT belong in dom{f}
    invalid = (X1_ == -1);
    f_(invalid) = NaN;

    % Plot f
    figure
    mesh(X1_,X2_,f_,'FaceColor','interp');
    hold on;
    [~,~,z0] = f(A,b,c,x_cvx);
    plot3(x_cvx(1), x_cvx(2), z0+1, 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 6);
    text(x_cvx(1) + 0.005, x_cvx(2) + 0.005, z0 + 1, 'x_*', 'FontSize', 12 ,'Color', 'k', 'FontWeight', 'bold');
    colorbar;
    xlabel('$x_1$','FontSize',14,'Interpreter','latex');
    ylabel('$x_2$','FontSize',14,'Interpreter','latex');
    hold off;

    % Plot f with visible barriers
    figure
    mesh(X1,X2,f_log_bar,'FaceColor','interp');
    hold on;
    [~,~,z0] = f(A,b,c,x_cvx);
    plot3(x_cvx(1), x_cvx(2), z0, 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 6);
    line ([x_cvx(1) x_cvx(1)], [x_cvx(2) x_cvx(2)], [z0 0], 'Color', 'r', 'LineWidth', 1.5, 'LineStyle', ':');
    text(x_cvx(1) + 0.005, x_cvx(2) + 0.005, z0 + 1, 'x_*', 'FontSize', 12 ,'Color', 'k', 'FontWeight', 'bold');
    colorbar;
    xlabel('$x_1$','FontSize',14,'Interpreter','latex');
    ylabel('$x_2$','FontSize',14,'Interpreter','latex');
    hold off;

    minf = min(f_(:));
    maxf = max(f_(:));

    ls = linspace(minf, maxf, 40);  % 40 evenly spaced lines

    % Plot the level sets of f
    figure
    contour(X1,X2,f_,ls);%,'ShowText','on');
    hold on;
    plot(x_cvx(1), x_cvx(2), 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 6);
    hold on;
    text(x_cvx(1) + 0.005, x_cvx(2) + 0.005, 'x_*', 'FontSize', 12 ,'Color', 'k', 'FontWeight', 'bold');
    xlabel('$x_1$','FontSize',14,'Interpreter','latex');
    ylabel('$x_2$','FontSize',14,'Interpreter','latex');
    % axis([min(0, x_cvx(1)) - 0.2 max(0, x_cvx(1)) + 0.2 min(0, x_cvx(2)) - 0.2 max(0, x_cvx(2)) + 0.2])
    hold off;
end

% 4. Perform GD w/BT-LS and x0 = 0
x0 = zeros(n,1);

alpha = 0.1;
beta = 0.7;
epsilon = 10^(-6);

fprintf('Now executing the gradient descent algorithm...\n');
[X_bt, F_bt] = Backtracking_LS_log_bar(A, b, c, x0, epsilon, alpha, beta, x_cvx);

if (round(F_bt(end),4) == round(cvx_optval,4))
    fprintf('GD with BT-LS found the same solution as cvx!\n\n');
end

fprintf('----------------\n');

% 5. Minimize f using the Newton algorithm
x0 = zeros(n,1);

alpha = 0.1;
beta = 0.7;
epsilon = 10^(-6);

fprintf('Now executing the Newton algorithm...\n');
[X_Nt, F_Nt] = Newton_log_bar(A, b, c, x0, epsilon, alpha, beta);

if (round(F_Nt(end),4) == round(cvx_optval,4))
    fprintf('Newton method found the same solution as cvx!\n\n');
end

fprintf('----------------\n');

% 6. Plot (f_GD - p_*) and (f_Nt - p_*) vs k
figure
semilogy([1:length(F_Nt)]', (F_Nt-F_Nt(end)), 'b', 'LineWidth', 1.5);
hold on;
semilogy([1:length(F_bt)]', (F_bt-F_bt(end)), 'r', 'LineWidth', 1.5);
grid on;
xlabel('$k$', 'FontSize', 13, 'Interpreter','latex');
ylabel('$log(f(\mathbf{x}_k)- p_*)$', 'FontSize', 13, 'Interpreter', 'latex');
legend('Newton', 'Backtracking LS');
axis tight