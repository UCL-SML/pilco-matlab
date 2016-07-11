%% calcCost.m
% *Summary:* Function to calculate the incurred cost and its standard deviation, 
% given a sequence of predicted state distributions and the cost struct
%
%   [L sL] = calcCost(cost, M, S)
%
% *Input arguments:*
%
%   cost               cost structure
%   M                  mean vectors of state trajectory (D-by-H matrix)
%   S                  covariance matrices at each time step (D-by-D-by-H)
%
% *Output arguments:*
%
%   L                  expected incurred cost of state trajectory
%   sL                 standard deviation of incurred cost
%
%
% Copyright (C) 2008-2013 by 
% Marc Deisenroth, Andrew McHutchon, Joe Hall, and Carl Edward Rasmussen.
%
% Last modified: 2013-01-23
%
%% High-Level Steps
% # Augment state distribution with trigonometric functions
% # Compute distribution of the control signal
% # Compute dynamics-GP prediction
% # Compute distribution of the next state
%

function [L sL] = calcCost(cost, M, S)
%% Code
H = size(M,2);                                             % horizon length
L = zeros(1,H); SL = zeros(1,H);

% for each time step, compute the expected cost and its variance
for h = 1:H
  [L(h),d1,d2,SL(h)]  = cost.fcn(cost, M(:,h), S(:,:,h));
end

sL = sqrt(SL);                                         % standard deviation

