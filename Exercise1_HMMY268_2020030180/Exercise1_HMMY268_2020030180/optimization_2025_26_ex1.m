% Convex Optimization 2025-26 %
% --- EXERCISE 1 --- %
% Nikolaos Malamas (2020030180)

clc; close all; clear all;

%% 
% 1.
step = 0.01;
x_max = 10;
x = [0:step:x_max];
f = 1./(x+1);

% (a) Done in paper - results in the report

% (b) We choose as x0 = [0,1,3] and plot the first and second order
% approximations at those points in common plots with the original f

% The closed forms of f_1 and f_2 have been derived by applying the
% given formulas in the exercise description. In the report, we
% we show in detail how the results were calculated.

% x0 = 0
f_1 = -x + 1;
f_2 = x.^2 - x + 1;

figure
plot(x,f);
hold on
plot(x,f_1);
hold on
plot(x,f_2);
x0 = 0;
y0 = f(x == x0);
xlabel('x','FontSize',12,'FontWeight','bold');
ylabel('y','FontSize',12,'FontWeight','bold');
plot(x0, 0, 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 6);
line([x0 x0], y0, 'Color', 'r', 'LineWidth', 2, 'LineStyle', ':');
text(x0 + 0.1, 3, 'x_0', 'FontSize', 12 ,'Color', 'k', 'FontWeight', 'bold');
grid on;
axis on
title(['$f(x)$ ' 'and its 1st and 2nd order Taylor approximations at ' '$x_0=0$'], 'Interpreter', 'latex');
legend('f(x)','f_1(x)','f_2(x)');

% x0 = 1
f_1 = -(1/4).*x + 3/4;
f_2 = (1/8).*x.^2 - (1/2).*x + 7/8;

figure
plot(x,f);
hold on
plot(x,f_1);
hold on
plot(x,f_2);
x0 = 1;
y0 = f(x == x0);
xlabel('x','FontSize',12,'FontWeight','bold');
ylabel('y','FontSize',12,'FontWeight','bold');
plot(x0, 0, 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 6);
line([x0 x0], [0 y0], 'Color', 'r', 'LineWidth', 1, 'LineStyle', ':');
text(x0 + 0.1, 0, 'x_0', 'FontSize', 12 ,'Color', 'k', 'FontWeight', 'bold');
grid on;
axis on
title(['$f(x)$ ' 'and its 1st and 2nd order Taylor approximations at ' '$x_0=1$'], 'Interpreter', 'latex');
legend('f(x)','f_1(x)','f_2(x)');

% x0 = 3
f_1 = -(1/16).*x + 7/16;
f_2 = (1/64).*x.^2 - (10/64).*x + 37/64;

figure
plot(x,f);
hold on
plot(x,f_1);
hold on
plot(x,f_2);
x0 = 3;
y0 = f(x == x0);
xlabel('x','FontSize',12,'FontWeight','bold');
ylabel('y','FontSize',12,'FontWeight','bold');
plot(x0, 0, 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 6);
line([x0 x0], [0 y0], 'Color', 'r', 'LineWidth', 1, 'LineStyle', ':');
text(x0 + 0.1, 0, 'x_0', 'FontSize', 12 ,'Color', 'k', 'FontWeight', 'bold');
grid on
axis on
title(['$f(x)$ ' 'and its 1st and 2nd order Taylor approximations at ' '$x_0=3$'], 'Interpreter', 'latex');
legend('f(x)','f_1(x)','f_2(x)');

%%
% 2.

% (a)
% Define a x* > 0 as a bound for the x1,x2 to be into [0,x*]
x_star = 10;
step = 0.1;

[x1,x2] = meshgrid(0:step:x_star);
f = 1./(x1+x2+1);

figure
mesh(x1,x2,f,"LineStyle","--","FaceColor","flat","EdgeColor","k","EdgeAlpha","0.2");
title('The plot of $f$ using the ${\tt mesh}$ function', 'Interpreter','latex');
xlabel('x_1','FontSize',12,'FontWeight','bold');
ylabel('x_2','FontSize',12,'FontWeight','bold');
zlabel('f(x_1,x_2)','FontSize',12,'FontWeight','bold');
colorbar

% (b)
no_levels = 15;
figure
contour(x1,x2,f,no_levels,LineWidth=2);
grid on;
axis([0 7 0 7]);
title('The plot of $f$ using the ${\tt contour}$ function','Interpreter','latex');
xlabel('x_1','FontSize',12,'FontWeight','bold');
ylabel('x_2','FontSize',12,'FontWeight','bold');
colorbar

% (c) Done in paper - results in the report
% The chosen x_0 = (0 0)
x01 = 0;
x02 = 0;
x0 = [x01;x02];

% (d)
% The 1st order Taylor approximation at x0
f1 = 1 - x1 - x2;

figure
mesh(x1,x2,f,"FaceColor","g","EdgeColor","k","EdgeAlpha","0.05",'LineWidth', 2.5);
hold on
mesh(x1,x2,f1,"FaceColor","b","EdgeColor","k","EdgeAlpha","0.05",'LineWidth', 2.5);
hold on
z0 = f(x1 == x01 & x2 == x02);
plot3(x01, x02, 0, 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 6);
line ([x01 0], [x02 0], [0 z0], 'Color', 'r', 'LineWidth', 1.5, 'LineStyle', ':');
text(x01 + 0.1, x02 + 0.1, 0, 'x_0', 'FontSize', 12 ,'Color', 'k', 'FontWeight', 'bold');
axis([0 2 0 2 -2 2]);
title('Common ${\tt mesh}$ of $f$ and its 1st order Taylor approximation $f_1$ at $x_0 = (0\ 0)$','Interpreter','latex');
legend('f(x)','f_1(x)');
xlabel('x_1','FontSize',12,'FontWeight','bold');
ylabel('x_2','FontSize',12,'FontWeight','bold');
zlabel('f(x_1,x_2)','FontSize',12,'FontWeight','bold');

% (e)
% The 2nd order Taylor approximation at x0
f2 = 1 - x1 - x2 + (x1+x2).^2;
figure
mesh(x1,x2,f,"FaceColor","g","EdgeColor","k","EdgeAlpha","0.1");
hold on
mesh(x1,x2,f2,"FaceColor","b","EdgeColor","k","EdgeAlpha","0.1");
hold on
z0 = f(x1 == x01 & x2 == x02);
plot3(x01, x02, 0, 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 6);
line([x01 0], [x02 0], [0 z0], 'Color', 'r', 'LineWidth', 1.5, 'LineStyle', ':');
text(x01 + 0.1, x02 + 0.1, 0, 'x_0', 'FontSize', 12 ,'Color', 'k', 'FontWeight', 'bold');
axis([0 2 0 2 -1 2]);
title('Common ${\tt mesh}$ of $f$ and its 2nd order Taylor approximation $f_2$ at $x_0 = (0\ 0)$','Interpreter','latex');
legend('f(x)','f_2(x)');
xlabel('x_1','FontSize',12,'FontWeight','bold');
ylabel('x_2','FontSize',12,'FontWeight','bold');
zlabel('f(x_1,x_2)','FontSize',12,'FontWeight','bold');

%%
% (5)
step = 0.01;
epsilon = 0.001;
xmax = 100;

% (a) Done in paper - results in the report

% (b) Done in paper - results in the report

% (c) Done in paper - results in the report

% (d) 
x = [epsilon:step:xmax];

% Define a color palette (10 distinct colors for 10 curves)
% using MATLAB's built-in color set
colors = turbo(10);   

figure
for i=1:5
    a = i;
    f_conv_gt1 = x.^a;
    f_concave = x.^(0.2*a);

    % Pick color indices for each pair of curves
    color_conv = colors(2*i-1, :);
    color_conc = colors(2*i, :);

    plot(x, f_conv_gt1, 'LineWidth', 2, 'Color', color_conv);
    hold on;
    plot(x, f_concave, 'LineWidth', 2, 'Color', color_conc);
    xlabel('x','FontSize',12,'FontWeight','bold');
    ylabel('f(x)','FontSize',12,'FontWeight','bold');
    title('Plot of $f(x)=x^a$ for $a \in[0,1]$ and $a\geq1$', 'Interpreter','latex');
    axis([0 3 0 5]);
    grid on;
    hold on;
end
legend({'f(x)=x^1','f(x)=x^{0.2}', ...
        'f(x)=x^2','f(x)=x^{0.4}', ...
        'f(x)=x^3','f(x)=x^{0.6}', ...
        'f(x)=x^4','f(x)=x^{0.8}', ...
        'f(x)=x^5','f(x)=x^{1.0}'}, ...
        'Location','northwest');

figure
for i=1:5
    a = i;
    f_conv_lt0 = x.^(-a);
    f_concave = x.^(0.2*a);

    % Pick color indices for each pair of curves
    color1 = colors(2*i-1, :);
    color2 = colors(2*i, :);

    plot(x, f_conv_lt0, 'LineWidth', 2, 'Color', color1);
    hold on;
    plot(x, f_concave, 'LineWidth', 2, 'Color', color2);
    xlabel('x','FontSize',12,'FontWeight','bold');
    ylabel('f(x)','FontSize',12,'FontWeight','bold');
    title('Plot of $f(x)=x^a$ for $a \in[0,1]$ and $a\leq0$', 'Interpreter','latex');
    axis([0 3 0 5]);
    grid on;
    hold on;
end
legend({'f(x)=x^{-1}','f(x)=x^{0.2}', ...
        'f(x)=x^{-2}','f(x)=x^{0.4}', ...
        'f(x)=x^{-3}','f(x)=x^{0.6}', ...
        'f(x)=x^{-4}','f(x)=x^{0.8}', ...
        'f(x)=x^{-5}','f(x)=x^{1.0}'}, ...
        'Location','northeast');

% (e)
xmax = 20;
[x1,x2]=meshgrid(-xmax:step:xmax);

f1 = sqrt(x1.^2 + x2.^2);
f2 = x1.^2 + x2.^2;

figure
mesh(x1, x2, f1, 'LineWidth', 2);
xlabel('x_1','FontSize',12,'FontWeight','bold');
ylabel('x_2','FontSize',12,'FontWeight','bold');
title('$f_1(\bf{x})=||\bf{x}||_2$,\ $n=2$', 'Interpreter','latex', 'FontSize', 14);
grid on;
axis on;

figure
mesh(x1, x2, f2, 'LineWidth', 2);
xlabel('x_1','FontSize',12,'FontWeight','bold');
ylabel('x_2','FontSize',12,'FontWeight','bold');
title('$f_2(\bf{x})=||\bf{x}||_2^2$,\ $n=2$', 'Interpreter','latex', 'FontSize', 14);
grid on;
axis on;


%%
% (6)
step = 0.1;

% Firstly, generate the positive definite matrix P
% as described in the exercise
A = randn(2,2);
P = A * A';
q = randn(2,1);
r = randn;

xmax = 20;

% Generate the mesh grid for x1 and x2
[X1, X2] = meshgrid(-xmax:step:xmax, -xmax:step:xmax);
x1 = [-xmax:step:xmax-step];
x2 = [-xmax:step:xmax-step];

f = (1/2).*(P(1,1)*X1.^2 + 2*P(1,2)*X1.*X2 + P(2,2)*X2.^2) ...
    + q(1)*X1 + q(2)*X2 + r;

x_star = -inv(P)*q;

figure
mesh(X1,X2,f);
hold on
z0 = (1/2).*(P(1,1)*x_star(1)^2 + 2*P(1,2)*x_star(1)*x_star(2) + P(2,2)*x_star(2)^2) ...
    + q(1)*x_star(1) + q(2)*x_star(2) + r;
plot3(x_star(1), x_star(2), 0, 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 8);
text(x_star(1) + 1, x_star(2) + 1, 150, 'x_*', 'FontSize', 14 ,'Color', 'k', 'FontWeight', 'bold');
xlabel('x_1','FontSize',12,'FontWeight','bold');
ylabel('x_2','FontSize',12,'FontWeight','bold');
title('The Quadratic Function $f(\mathbf{x})=\frac{1}{2}\mathbf{x}^T\mathbf{P}\mathbf{x}+\mathbf{q}^T\mathbf{x}+r$,\ in {\tt mesh} plot','Interpreter','latex');
grid on;
colorbar

figure
contour(X1, X2, f, no_levels, LineWidth=2);
hold on;
plot(x_star(1), x_star(2), 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 6);
text(x_star(1) + 0.5, x_star(2) + 0.5, 'x_*', 'FontSize', 12 ,'Color', 'k', 'FontWeight', 'bold');
xlabel('x_1','FontSize',12,'FontWeight','bold');
ylabel('x_2','FontSize',12,'FontWeight','bold');
title('The Quadratic Function $f(\mathbf{x})=\frac{1}{2}\mathbf{x}^T\mathbf{P}\mathbf{x}+\mathbf{q}^T\mathbf{x}+r$,\ in {\tt contour} plot','Interpreter','latex');
grid on;
colorbar

