%% learnPolicy.m
% *Summary:* Script to perform the policy search
%
% Copyright (C) 2008-2013 by
% Marc Deisenroth, Andrew McHutchon, Joe Hall, and Carl Edward Rasmussen.
%
% Last modified: 2013-03-06
%
%% High-Level Steps
% # Learn the policy (call optimizer)
% # Predict trajectory from a sampled start state and compute expected cost

%% Code

% 1. Update the policy
opt.fh = 1;
[policy.p fX3] = minimize(policy.p, 'value', opt, mu0Sim, S0Sim, ...
  dynmodel, policy, plant, cost, H);

% (optional) Plot overall optimization progress
if exist('plotting', 'var') && isfield(plotting, 'verbosity') ...
    && plotting.verbosity > 1
  if ~ishandle(2); figure(2); else set(0,'CurrentFigure',2); end
  hold on; plot(fX3); drawnow; 
  xlabel('line search iteration'); ylabel('function value')
end

% 2. Predict state trajectory from p(x0) and compute cost trajectory
[M{j} Sigma{j}] = pred(policy, plant, dynmodel, mu0Sim(:,1), S0Sim, H);
[fantasy.mean{j} fantasy.std{j}] = ...
  calcCost(cost, M{j}, Sigma{j}); % predict cost trajectory

% (optional) Plot predicted immediate costs (as a function of the time steps)
if exist('plotting', 'var') && isfield(plotting, 'verbosity') ...
    && plotting.verbosity > 0
  if ~ishandle(3); figure(3); else set(0,'CurrentFigure',3); end
  clf(3); errorbar(0:H,fantasy.mean{j},2*fantasy.std{j}); drawnow;
  xlabel('time step'); ylabel('immediate cost');
end