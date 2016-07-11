%% getPlotDistr_pendubot.m
% *Summary:* Compute means and covariances of the Cartesian coordinates of
% the tips both the inner and outer pendulum assuming that the joint state
% $x$ of the cart-double-pendulum system is Gaussian, i.e., $x\sim N(m, s)$
%
%
%     function [M1, S1, M2, S2] = getPlotDistr_pendubot(m, s, ell1, ell2)
%
%
%
% *Input arguments:*
%
%   m       mean of full state                                    [6 x 1]
%   s       covariance of full state                              [6 x 6]
%   ell1    length of inner pendulum
%   ell2    length of outer pendulum
%
%   Note: this code assumes that the following order of the state:
%          1: pend1 angular velocity,
%          2: pend2 angular velocity, 
%          3: pend1 angle, 
%          4: pend2 angle
%
% *Output arguments:*
%
%   M1      mean of tip of inner pendulum                         [2 x 1]
%   S1      covariance of tip of inner pendulum                   [2 x 2]
%   M2      mean of tip of outer pendulum                         [2 x 1]
%   S2      covariance of tip of outer pendulum                   [2 x 2]
%
% Copyright (C) 2008-2013 by 
% Marc Deisenroth, Andrew McHutchon, Joe Hall, and Carl Edward Rasmussen.
%
% Last modification: 2013-03-27
%
%% High-Level Steps
% # Augment input distribution to complex angle representation
% # Compute means of tips of pendulums (in Cartesian coordinates)
% # Compute covariances of tips of pendulums (in Cartesian coordinates)

function [M1, S1, M2, S2] = getPlotDistr_pendubot(m, s, ell1, ell2)
%% Code

% 1. Augment input distribution
[m1 s1 c1] = gTrig(m, s, [3 4], [ell1, ell2]); % map input through sin/cos
m1 = [m; m1];        % mean of joint
c1 = s*c1;           % cross-covariance between input and prediction
s1 = [s c1; c1' s1]; % covariance of joint

% 2. Compute means of tips of pendulums (in Cartesian coordinates)
M1 = [-m1(5); m1(6)];                 % [-l*sin(t1), l*cos(t1)]
M2 = [-m1(5) + m1(7); m1(6) + m1(8)]; % [-l*(sin(t1)-sin(t2)),l*(cos(t1)+cos(t2))]

% 2. Put covariance matrices together (Cart. coord.)
% first set of coordinates (tip of 1st pendulum)
s11 = s1(5,5);
s12 = -s1(5,6);
s22 = s1(6,6);
S1 = [s11 s12; s12' s22];

% second set of coordinates (tip of 2nd pendulum)
s11 = s1(5,5) + s1(7,7) - s1(5,7) - s1(7,5);    % ell1*sin(t1) + ell2*sin(t2)
s22 = s1(6,6) + s1(8,8) + s1(6,8) + s1(8,6);    % ell1*cos(t1) + ell2*cos(t2)
s12 = -(s1(5,6) + s1(5,8) + s1(7,6) + s1(7,8)); 
S2 = [s11 s12; s12' s22];

% make sure we have proper covariances (sometimes numerical problems occur)
try
  chol(S1); 
catch
  warning('matrix S1 not pos.def. (getPlotDistr)');
  S1 = S1 + (1e-6 - min(eig(S1)))*eye(2);
end

try
  chol(S2); 
catch
  warning('matrix S2 not pos.def. (getPlotDistr)');
  S2 = S2 + (1e-6 - min(eig(S2)))*eye(2);
end