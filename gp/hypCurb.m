%% hypCurb.m
% *Summary:* Wrapper for GP training (via gpr.m), penalizing large SNR and 
% extreme length-scales to avoid numerical instabilities
%
%     function [f df] = hypCurb(lh, covfunc, x, y, curb)
%
% *Input arguments:*
%
%   lh       log-hyper-parameters                                   [D+2 x  E ]
%   covfunc  covariance function, e.g., 
%                               covfunc = {'covSum', {'covSEard', 'covNoise'}};
%   x        training inputs                                        [ n  x  D ]
%   y        training targets                                      [ n  x  E ]
%   curb     (optional) parameters to penalize extreme hyper-parameters
%     .ls    length-scales
%     .snr   signal-to-noise ratio (try to keep it below 500)
%     .std   additional parameter required for length-scale penalty 
%
% *Output arguments:*
%
%   f        penalized negative log-marginal likelihood
%   df       derivative of penalized negative log-marginal likelihood wrt
%            GP log-hyper-parameters
%
%
% Copyright (C) 2008-2013 by
% Marc Deisenroth, Andrew McHutchon, Joe Hall, and Carl Edward Rasmussen.
%
% Last modified: 2011-12-19
%
%% High-Level Steps
% # Compute the negative log-marginal likelihood (plus derivatives)
% # Add penalties and change derivatives accordingly

function [f df] = hypCurb(lh, covfunc, x, y, curb)
%% Code
if nargin < 5, curb.snr = 500; curb.ls = 100; curb.std = 1; end   % set default

p = 30;                                                         % penalty power
D = size(x,2);
if size(lh,1) == 3*D+2; li = 1:2*D; sfi = 2*D+1:3*D+1; % 1D and DD terms
elseif size(lh,1) == 2*D+1; li = 1:D; sfi = D+1:2*D;   % Just 1D terms
elseif size(lh,1) == D+2; li = 1:D; sfi = D+1;         % Just DD terms
else error('Incorrect number of hyperparameters'); 
end

ll = lh(li); lsf = lh(sfi); lsn = lh(end);

% 1) compute the negative log-marginal likelihood (plus derivatives)
[f df] = gpr(lh, covfunc, x, y);                              % first, call gpr

% 2) add penalties and change derivatives accordingly
f = f + sum(((ll - log(curb.std'))./log(curb.ls)).^p);   % length-scales
df(li) = df(li) + p*(ll - log(curb.std')).^(p-1)/log(curb.ls)^p;

f = f + sum(((lsf - lsn)/log(curb.snr)).^p); % signal to noise ratio
df(sfi) = df(sfi) + p*(lsf - lsn).^(p-1)/log(curb.snr)^p;
df(end) = df(end) - p*sum((lsf - lsn).^(p-1)/log(curb.snr)^p);
