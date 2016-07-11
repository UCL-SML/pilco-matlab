% plot the squashing function

clear all; close all;
header;


x = linspace(-5/2*pi,5/2*pi,1000);

squash = @(x) (9*sin(x) + sin(3.*x))./8;


y = squash(x);

% plot a few periods 
figure(1);
clf; hold on;
plot(x,y);
xlabel('$x$', 'interpreter', 'latex');
ylabel('$\sigma(x)$', 'interpreter', 'latex');
axis tight;

set(gcf,'PaperSize', [12 6]);
set(gcf,'PaperPosition',[0.1 0.1 12 6]);
print_pdf('../figures/squashing_fct');

% plot a single period
axis([-pi/2, pi/2, -1, 1])
print_pdf('../figures/squashing_fct2');