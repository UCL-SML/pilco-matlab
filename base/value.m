%% value.m
% *Summary:* Compute expected (discounted) cumulative cost for a given (set of) initial
% state distributions
%
%     function [J, dJdp] = value(p, m0, S0, dynmodel, policy, plant, cost, H)
%
% *Input arguments:*
%
%   p            policy parameters chosen by minimize
%   policy       policy structure
%     .fcn       function which implements the policy
%     .p         parameters passed to the policy
%   m0           matrix (D by k) of initial state means
%   S0           covariance matrix (D by D) for initial state
%   dynmodel     dynamics model structure
%   plant        plant structure
%   cost         cost function structure
%     .fcn       function handle to the cost
%     .gamma     discount factor
%  H             length of prediction horizon
%
% *Output arguments:*
%
%  J             expected cumulative (discounted) cost
%  dJdp          (optional) derivative of J wrt the policy parameters
%
% Copyright (C) 2008-2013 by 
% Marc Deisenroth, Andrew McHutchon, Joe Hall, and Carl Edward Rasmussen.
%
% Last modification: 2013-03-21
%
%% High-Level Steps
% # Compute distribution of next state
% # Compute corresponding expected immediate cost (discounted)
% # At end of prediction horizon: sum all immediate costs up

function [J, dJdp] = value(p, m0, S0, dynmodel, policy, plant, cost, H)
%% Code

policy.p = p;            % overwrite policy.p with new parameters from minimize
p = unwrap(policy.p); dp = 0*p;
m = m0; S = S0; L = zeros(1,H);

if nargout <= 1                                       % no derivatives required
  
  for t = 1:H                                  % for all time steps in horizon
    [m, S] = plant.prop(m, S, plant, dynmodel, policy);      % get next state
    L(t) = cost.gamma^t.*cost.fcn(cost, m, S);     % expected discounted cost
  end
  
else                                               % otherwise, get derivatives
  
  dmOdp = zeros([size(m0,1), length(p)]);
  dSOdp = zeros([size(m0,1)*size(m0,1), length(p)]);
  
  for t = 1:H                                  % for all time steps in horizon
    [m, S, dmdmO, dSdmO, dmdSO, dSdSO, dmdp, dSdp] = ...
      plant.prop(m, S, plant, dynmodel, policy); % get next state
    
    dmdp = dmdmO*dmOdp + dmdSO*dSOdp + dmdp;
    dSdp = dSdmO*dmOdp + dSdSO*dSOdp + dSdp;
    
    [L(t), dLdm, dLdS] = cost.fcn(cost, m, S);              % predictive cost
    L(t) = cost.gamma^t*L(t);                                      % discount
    dp = dp + cost.gamma^t*( dLdm(:)'*dmdp + dLdS(:)'*dSdp )';
    
    dmOdp = dmdp; dSOdp = dSdp;                                 % bookkeeping
  end
  
end

J = sum(L); dJdp = rewrap(policy.p, dp);
