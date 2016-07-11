%% gaussian.m
% *Summary:* Generate n samples from a Gaussian $p(x)=\mathcal N(m,S).
% Sampling is based on the Cholesky factorization of the covariance matrix S
%
%    function x = gaussian(m, S, n)
%
% *Input arguments:*
%
%   m      mean of Gaussian                                         [D x 1]
%   S      covariance matrix of Gaussian                            [D x D]
%   n      (optional) number of samples; default: n=1
%
% *Output arguments:*
%
%   x      matrix of samples from Gaussian                          [D x n]
%
%
% Copyright (C) 2008-2013 by 
% Marc Deisenroth, Andrew McHutchon, Joe Hall, and Carl Edward Rasmussen.
%
% Last modified: 2013-03-21

function x = gaussian(m, S, n)
%% Code

if nargin < 3; n = 1; end

x = bsxfun(@plus, m(:), chol(S)'*randn(size(S,2),n));
