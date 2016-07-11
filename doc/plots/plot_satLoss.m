% plot the saturating cost function

clear all; close all;
header;

addpath ('../../loss/');
addpath ('../../util/');

cost.fcn = @lossSat;
cost.z = 0;              % target state
cost.width = 0.5;


cost.W = 1/(cost.width).^2; % weight matrix

x = linspace(0,3,200)';  % array of distances

c = zeros(1,length(x));

% evaluate cost for all distances x
for i = 1:length(x)
  c(i) = cost.fcn(cost, x(i), 0);
end

% evaluate cost at 2-sigma bound
c2 = cost.fcn(cost, 2*cost.width, 0);

figure; 
hold on
plot(x,c); % plot cost
plot(2*cost.width, c2, 'ro', 'markersize', 8); % plot 2-sigma cost value
plot(linspace(0,2*cost.width,100), c2*ones(1,100), 'k-'); % plot line
text(cost.width/4, c2+0.05, '$2\texttt{cost.width}$','interpreter', 'latex')
xlabel('Distance to target')
ylabel('Cost')
axis([x(1) x(end) 0 1.01]);

set(gcf,'PaperSize', [12 6]);
set(gcf,'PaperPosition',[0.1 0.1 12 6]);
print_pdf('../figures/satLossPlot');
