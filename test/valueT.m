%% valueT.m
% *Summary:* Test derivatives of the propagate function, which computes the 
% mean and the variance of the successor state distribution, assuming that the
% current state is Gaussian distributed with mean m and covariance matrix
% s.
%
%   [d dy dh] = valueT(p, delta, m, s, dynmodel, policy, plant, cost, H)
%
%
% *Input arguments:*
%
%   p          policy parameters (can be a structure)
%     .<>      fields that contain the policy parameters (nothing else)
%   m          mean of the input distribution
%   s          covariance of the input distribution
%   dynmodel   GP dynamics model (structure)
%   policy     policy structure
%   plant      plant structure
%   cost       cost structure
%   H          prediction horizon
%   delta      (optional) finite difference parameter. Default: 1e-4
%
%
% *Output arguments:*
%
%   dd         relative error of analytical vs. finite difference gradient
%   dy         analytical gradient
%   dh         finite difference gradient
%
% Copyright (C) 2008-2013 by 
% Marc Deisenroth, Andrew McHutchon, Joe Hall, and Carl Edward Rasmussen.
%
% Last modified: 2013-03-21

function [d dy dh] = valueT(p, m, s, dynmodel, policy, plant, cost, H, delta)
%% Code 

if nargin < 9; delta = 1e-4; end
if nargin < 8; H = 4; end

% call checkgrad directly
[d dy dh] = checkgrad('value',p,delta,m,s,dynmodel,policy,plant,cost,H);
