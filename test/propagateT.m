%% propagateT.m
% *Summary:* Test derivatives of the propagate function, which computes the 
% mean and the variance of the successor state distribution, assuming that the
% current state is Gaussian distributed with mean m and covariance matrix
% s.
%
%   [dd dy dh] = propagateT(deriv, plant, dynmodel, policy, m, s, delta)
%
%
% *Input arguments:*
%
%   deriv    desired derivative. options:
%        (i)    'dMdm' - derivative of the mean of the predicted state
%                wrt the mean of the input distribution
%        (ii)   'dMds' - derivative of the mean of the predicted state 
%                wrt the variance of the input distribution
%        (iii)  'dMdp' - derivative of the mean of the predicted state 
%                wrt the controller parameters
%        (iv)   'dSdm' - derivative of the variance of the predicted state 
%                wrt the mean of the input distribution
%        (v)    'dSds' - derivative of the variance of the predicted state 
%                wrt the variance of the input distribution
%        (vi)   'dSdp' - derivative of the variance of the predicted state 
%                wrt the controller parameters
%        (vii)  'dCdm' - derivative of the inv(s)*(covariance of the input and 
%                the predicted state) wrt the mean of the input distribution
%        (viii) 'dCds' - derivative of the inv(s)*(covariance of the input and 
%                the predicted state) wrt the variance of the input distribution
%        (ix)   'dCdp' - derivative of the inv(s)*(covariance of the input and 
%                the predicted state) wrt the controller parameters
%   plant      plant structure
%     .poli    state indices: policy inputs
%     .dyno    state indices: predicted variables
%     .dyni    state indices: inputs to ODE solver
%     .difi    state indices that are learned via differences
%     .angi    state indices: angles
%     .poli    state indices: policy inputs
%     .prop    function handle to function responsible for state
%              propagation. Here: @propagated
%   dynmodel   GP dynamics model (structure)
%     .hyp     log-hyper parameters
%     .inputs  training inputs
%     .targets training targets
%   policy     policy structure
%     .maxU    maximum amplitude of control
%     .fcn     function handle to policy 
%     .p       struct of policy parameters
%     .p.<>    policy-specific parameters are stored here
%   m          mean of the input distribution
%   s          covariance of the input distribution
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


function [dd dy dh] = propagateT(deriv, plant, dynmodel, policy, m, s, delta)
%% Code
if nargin < 2,
  randn('seed',24)
  nn = 10;
  E = 4; F = 3;
  m = randn(D,1);
  s = randn(D); s = s*s';
  
  % Plant ----------------------------------------------------------------
  plant.poli = 1:E;
  plant.dyno = 1:E;
  plant.dyni = 1:E;
  plant.difi = 1:E;
  plant.angi = [];
  plant.prop = @propagated;
  
  % Policy ---------------------------------------------------------------
  policy.p.w = randn(F,E);
  policy.p.b = randn(F,1);
  policy.fcn = @(policy,m,s)conCat(@conlin,@gSat,policy,m,s);
  policy.maxU = 20*ones(1,F);
  
  % Dynamics -------------------------------------------------------------
  dynmodel.hyp = zeros(1,E);
  dynmodel.inputs  = randn(nn,E+F);
  dynmodel.targets = randn(nn,E);
  
end

if nargin < 7; delta = 1e-4; end


D = length(m);

switch deriv
  
  case 'dMdm'      
      [dd dy dh] = checkgrad(@propagateT0, m, delta, plant, dynmodel, policy, s);

  case 'dSdm'
      [dd dy dh] = checkgrad(@propagateT1, m, delta, plant, dynmodel, policy, s);

  case 'dMds'
      [dd dy dh] = checkgrad(@propagateT2, s(tril(ones(D))==1), delta, plant, ...
                                                        dynmodel, policy, m);
   
  case 'dSds'
      [dd dy dh] = checkgrad(@propagateT3, s(tril(ones(D))==1), delta, plant, ...
                                                        dynmodel, policy, m);
 
  case 'dMdp'
      p = unwrap(policy.p);
      [dd dy dh] = checkgrad(@propagateT4, p, delta, plant, dynmodel, policy, m, s);
 
  case 'dSdp'
      p = unwrap(policy.p);
      [dd dy dh] = checkgrad(@propagateT5, p, delta, plant, dynmodel, policy, m, s);
end


function [f, df] = propagateT0(m, plant, dynmodel, policy, s)       % dMdm
if nargout == 1
  M = plant.prop(m, s, plant, dynmodel, policy);
else
  [M, S, dMdm] = plant.prop(m, s, plant, dynmodel, policy);
  df = permute(dMdm,[1,3,2]);
end
f = M;

function [f, df] = propagateT1(m, plant, dynmodel, policy, s)       % dSdm
if nargout == 1
  [M, S] = plant.prop(m, s, plant, dynmodel, policy);
else
  [M, S, dMdm, dSdm] = plant.prop(m, s, plant, dynmodel, policy);
  dd = length(M); df = reshape(dSdm,dd,dd,[]);
end
f = S;

function [f, df] = propagateT2(s, plant, dynmodel, policy, m)       % dMds
d = length(m);
ss(tril(ones(d))==1) = s; ss = reshape(ss,d,d); ss = ss + ss' - diag(diag(ss));
if nargout == 1
  M = plant.prop(m, ss, plant, dynmodel, policy);
else
  [M, S, dMdm, dSdm, dMds] = plant.prop(m, ss, plant, dynmodel, policy);
  dd = length(M); dMds = reshape(dMds,dd,d,d); df = zeros(dd,1,d*(d+1)/2);
  for i=1:dd; 
        dMdsi(:,:) = dMds(i,:,:); dMdsi = dMdsi + dMdsi'-diag(diag(dMdsi)); 
        df(i,:) = dMdsi(tril(ones(d))==1);
  end
end
f = M;

function [f, df] = propagateT3(s, plant, dynmodel, policy, m)       % dSds
d = length(m);
ss(tril(ones(d))==1) = s; ss = reshape(ss,d,d); ss = ss + ss' - diag(diag(ss));
if nargout == 1
  [M, S] = plant.prop(m, ss, plant, dynmodel, policy);
else
  [M, S, dMdm, dSdm, dMds, dSds] = ...
                               plant.prop(m, ss, plant, dynmodel, policy);
  dd = length(M); dSds = reshape(dSds,dd,dd,d,d); df = zeros(dd,dd,d*(d+1)/2);
    for i=1:dd; for j=1:dd
        dSdsi(:,:) = dSds(i,j,:,:); dSdsi = dSdsi+dSdsi'-diag(diag(dSdsi)); 
        df(i,j,:) = dSdsi(tril(ones(d))==1);
    end; end
end
f = S;

function [f, df] = propagateT4(p, plant, dynmodel, policy, m, s)    % dMdp
policy.p = rewrap(policy.p, p);
if nargout == 1
  M = plant.prop(m, s, plant, dynmodel, policy);
else
  [M, S, dMdm, dSdm, dMds, dSds, dMdp] = plant.prop(m, s, plant, dynmodel, policy);
  df = dMdp;
end
f = M;

function [f, df] = propagateT5(p, plant, dynmodel, policy, m, s)    % dSdp
policy.p = rewrap(policy.p, p);
if nargout == 1
  [M, S] = plant.prop(m, s, plant, dynmodel, policy);
else
  [M, S, dMdm, dSdm, dMds, dSds, dMdp, dSdp] = ...
                                plant.prop(m, s, plant, dynmodel, policy);
  df = dSdp;
end
f = S;
