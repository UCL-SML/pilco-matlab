header;
clear all; close all
addpath(genpath('../../'));
randn('state',7); rand('state', 1);

covfunc = {'covSum', {'covSEard', 'covNoise'}};        % specify ARD covariance


n = 50;

Axis = [-5 5 -1.5 2];

dynmodel.inputs = rand(n,1)*10 - 5;
dynmodel.hyp = [0.2, 1, 0.01]';
dynmodel.target = 1.5*randn(n,1);

xx = linspace(-5,5,101)';
xx2 = linspace(-pi/2,pi/2, 101);

prelPolicy = gpr(dynmodel.hyp, covfunc, dynmodel.inputs, dynmodel.target, xx);

figure(1); clf; hold on;
plot(xx, prelPolicy);
plot(xx, ones(size(xx)), 'k--');
plot(xx, -ones(size(xx)), 'k--');

xlabel('$x$', 'interpreter', 'latex')
ylabel('$\tilde\pi(x)$', 'interpreter', 'latex')
axis(Axis)
set(gcf,'PaperSize', [12 6]); 
set(gcf,'PaperPosition',[0.1 0.1 12 6]);
print_pdf('../figures/preliminary_policy');

for i = 1:size(xx,1)
  squashedPolicy(i) = gSat(prelPolicy(i), 0, 1);
  sqashingFct(i) = gSat(xx2(i), 0, 1);
end

figure(2); clf; hold on;
plot(xx, squashedPolicy);
plot(xx, ones(size(xx)), 'k--');
plot(xx, -ones(size(xx)), 'k--');
xlabel('$x$', 'interpreter', 'latex')
ylabel('$\pi(x)$','interpreter', 'latex')
axis(Axis)
set(gcf,'PaperSize', [12 6]); set(gcf,'PaperPosition',[0.1 0.1 12 6]);
print_pdf('../figures/squashed_preliminary_policy');
