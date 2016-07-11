%% dynamics_cp.m
% *Summary:* Implements ths ODE for simulating the cart-pole dynamics.
%
%    function dz = dynamics_cp(t, z, f)
%
%
% *Input arguments:*
%
%		t     current time step (called from ODE solver)
%   z     state                                                    [4 x 1]
%   f     (optional): force f(t)
%
% *Output arguments:*
%
%   dz    if 3 input arguments:      state derivative wrt time
%         if only 2 input arguments: total mechanical energy
%
%
% Note: It is assumed that the state variables are of the following order:
%       x:        [m]     position of cart
%       dx:       [m/s]   velocity of cart
%       dtheta:   [rad/s] angular velocity
%       theta:    [rad]   angle
%
%
% A detailed derivation of the dynamics can be found in:
%
% M.P. Deisenroth:
% Efficient Reinforcement Learning Using Gaussian Processes, Appendix C,
% KIT Scientific Publishing, 2010.
%
% Copyright (C) 2008-2013 by
% Marc Deisenroth, Andrew McHutchon, Joe Hall, and Carl Edward Rasmussen.
%
% Last modified: 2013-03-08

function dz = dynamics_cp(t,z,f)
%% Code

l = 0.5;  % [m]      length of pendulum
m = 0.5;  % [kg]     mass of pendulum
M = 0.5;  % [kg]     mass of cart
b = 0.1;  % [N/m/s]  coefficient of friction between cart and ground
g = 9.82; % [m/s^2]  acceleration of gravity

if nargin==3
  dz = zeros(4,1);
  dz(1) = z(2);
  dz(2) = ( 2*m*l*z(3)^2*sin(z(4)) + 3*m*g*sin(z(4))*cos(z(4)) ...
          + 4*f(t) - 4*b*z(2) )/( 4*(M+m)-3*m*cos(z(4))^2 );
  dz(3) = (-3*m*l*z(3)^2*sin(z(4))*cos(z(4)) - 6*(M+m)*g*sin(z(4)) ...
          - 6*(f(t)-b*z(2))*cos(z(4)) )/( 4*l*(m+M)-3*m*l*cos(z(4))^2 );
  dz(4) = z(3);
else
  dz = (M+m)*z(2)^2/2 + 1/6*m*l^2*z(3)^2 + m*l*(z(2)*z(3)-g)*cos(z(4))/2;
end