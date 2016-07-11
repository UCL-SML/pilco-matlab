%% lossLin.m
% *Summary:* Function to compute the expected loss and its derivatives, given an 
% input distribution, under a linear loss function: L = a^T(x - b). Note, this
% loss function can return negative loss.
%
%   [L dLdm dLds S dSdm dSds C dCdm dCds] = lossLin(cost,m,s)
%
% *Input arguments:*
%
%    cost
%      .a      gradient of linear loss function,                       [D x 1]
%      .b      targets, the value of x for which there is zero loss    [D x 1]
%    m         mean of input distribution                              [D x 1]
%    s         covariance matrix of input distribution                 [D x D]
%
% *Output arguments:*
%
%   L               expected loss                                  [1   x    1 ]
%   dLdm            derivative of L wrt input mean                 [1   x    D ]
%   dLds            derivative of L wrt input covariance           [1   x   D^2] 
%   S               variance of loss                               [1   x    1 ]
%   dSdm            derivative of S wrt input mean                 [1   x    D ]
%   dSds            derivative of S wrt input covariance           [1   x   D^2]    
%   C               inv(S) times input-output covariance           [D   x    1 ]   
%   dCdm            derivative of C wrt input mean                 [D   x    D ]  
%   dCds            derivative of C wrt input covariance           [D   x   D^2]  
%
% Copyright (C) 2008-2013 by
% Marc Deisenroth, Andrew McHutchon, Joe Hall, and Carl Edward Rasmussen.
%
% Last modified: 2013-03-05
%
%% High-Level Steps
% # Expected cost
% # Variance of cost
% # inv(s)* cov(x,L)

function [L dLdm dLds S dSdm dSds C dCdm dCds] = lossLin(cost,m,s)
%% Code

a = cost.a(:); b = cost.b(:); D = length(m);
if length(a) ~= D || length(b) ~= D;
    error('a or b not the same length as m'); end

% 1) Mean
L = a'*(m - b);
dLdm = a';
dLds = zeros(D);

% 2) Variance
S = a'*s*a;
dSdm = zeros(1,D);
dSds = a*a';

% 3) inv(s) * input-output covariance cov(x,L)
C = a;
dCdm = zeros(D); dCds = zeros(D,D^2);

