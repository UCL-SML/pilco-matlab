clear all; close all;
try
  rd = '../../../';
  addpath([rd 'base'],[rd 'util'],[rd 'gp'],[rd 'control'],[rd 'loss']);
  addpath ../
catch
end

iter = ['1', '2', '3', '4', '5'];
success = zeros(12,length(iter));

%HH = 100;

for horizon = 1:20
    for numTrials = 1:length(iter)
        filename = ['doublepend', iter(numTrials), '_', num2str(horizon), '_H40.mat'];
        load(filename)
        % regenerate a rollout with a longer horizon
        %xx = rollout(gaussian(mu0, S0), policy, HH, plant, cost);
        %disp(xx)
        fprintf('Iteration: %i', horizon);
        xy = mod(xx(end-5:end,3:4), 2*pi);
        if min(min( xy  < deg2rad(10) | xy > deg2rad(350)))
            success(horizon,numTrials) = 1;
            fprintf 'Success\n'
        end
    end
end


