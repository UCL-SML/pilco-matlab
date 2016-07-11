%% dynamics_unicycle.m
% *Summary:* Implements ths ODE for simulating the cart-pole dynamics. 
%
%    function dz = dz = dynamics_unicycle(t, z, V, U)
%
%
% *Input arguments:*
%
%		t     current time step (called from ODE solver)
%   z     state                                                    [12 x 1]
%   V     torque applied to the flywheel
%   U     torque applied to the wheel
%
% *Output arguments:*
%   
%   dz    state derivative wrt time
%
%
% Note: It is assumed that the state variables are of the following order:
% state: z = [dtheta, dphi, dpsiw, dpsif, dspit,
%             x,  y,  theta,  phi,  psiw,  psif,  psit]
%
%   theta: tilt of the unicycle
%   phi: orientation of the unicycle
%   psiw: angle of wheel (rotation)
%   psif: angle of fork
%   psit: angle of turntable (rotation)
%
%       dtheta   angular velocity of tilt of the unicycle
%       dphi     angular velocity of orientation of the unicycle
%       dpsiw    angular velocity of wheel
%       dpsif    angular velocity of fork
%       dpsit    angular velocity of turntable
%       x        x-position of contact point in plane
%       y        y-position of contact point in plane
%       theta    tilt of the unicycle
%       phi      orientation of the unicycle
%       psiw     angle of wheel (rotation)
%       psif     angle of fork
%       psit     angle of turntable (rotation)    
%
%
% Copyright (C) 2008-2013 by
% Marc Deisenroth, Andrew McHutchon, Joe Hall, and Carl Edward Rasmussen,
% based on derivations by David Forster
%
% Last modified: 2013-03-18

function dz = dynamics_unicycle(t, z, V, U)
%% Code

T = 0; % no friction

dtheta = z(1); dphi = z(2); dpsiw = z(3); dpsif = z(4); dpsit = z(5); x = z(6);
y = z(7); theta = z(8); phi = z(9); psiw = z(10); psif = z(11); psit = z(12);
 
clear psit psiw   % dynamics can't possibly depend on these
 
% plant characteristics 

mt = 10.0;   % turntable mass
mw =  1.0;   % wheel mass
mf = 23.5;   % frame mass
rw =  0.225; % wheel radius 
rf =  0.54;  % frame center of mass to wheel
rt =  0.27;  % frame centre of mass to turntable
r =  rf+rt;  % distance wheel to turntable
Cw = 0.0484; % moment of inertia of wheel around axle
Aw = 0.0242; % moment of inertia of wheel perpendicular to axle
Cf = 0.8292; % moment of inertia of frame
Bf = 0.4608; % moment of inertia of frame
Af = 0.4248; % moment of inertia of frame
Ct = 0.2;    % moment of inertia of turntable around axle
At = 1.3;    % moment of inertia of turntable perpendicular to axle
g = 9.82;    % acceleration of gravity

st = sin(theta); ct = cos(theta); sf = sin(psif); cf = cos(psif);

A = [                                                                                                                     -Ct*sf                                                                                                                                                                Ct*cf*ct                                                                                                  0                                                 0    Ct;
                                                                                                                               0                                                                                          Cw*st+At*st-rf*(-mf*(st*rf+cf*st*rw)-mt*(st*r+cf*st*rw))+rt*mt*(st*r+cf*st*rw)                                                                          -cf*rw*(rf*(mf+mt)+rt*mt)                    -Cw-At-rf*(mf*rf+mt*r)-rt*mt*r     0;
                                            cf*(-Af*sf-Ct*sf)-sf*(-Bf*cf-At*cf+rf*(-mf*(cf*rf+rw)-mt*(cf*r+rw))-rt*mt*(cf*r+rw))                                                                         Aw*ct+cf*(Af*cf*ct+Ct*cf*ct)-sf*(-Bf*sf*ct-At*sf*ct+rf*(-mf*sf*ct*rf-mt*sf*ct*r)-rt*mt*sf*ct*r)                                                                                                  0                                                 0 Ct*cf;
  -Aw-rw*(mf*(cf*rf+rw)+mw*rw+mt*(cf*r+rw))+sf*(-Af*sf-Ct*sf)+cf*(-Bf*cf-At*cf+rf*(-mf*(cf*rf+rw)-mt*(cf*r+rw))-rt*mt*(cf*r+rw))                                                  -rw*(mt*sf*ct*r+mf*sf*ct*rf)+sf*(Af*cf*ct+Ct*cf*ct)+cf*(-Bf*sf*ct-At*sf*ct+rf*(-mf*sf*ct*rf-mt*sf*ct*r)-rt*mt*sf*ct*r)                                                                                                  0                                                 0 Ct*sf;
                                                                                                                               0 2*Cw*st+At*st-rf*(-mt*(st*r+cf*st*rw)-mf*(st*rf+cf*st*rw))+rt*mt*(st*r+cf*st*rw)+rw*(mw*st*rw+sf*(mf*sf*st*rw+mt*sf*st*rw)+cf*(mt*(st*r+cf*st*rw)+mf*(st*rf+cf*st*rw))) -Cw-rt*mt*cf*rw+rw*(-mw*rw+sf*(-mf*sf*rw-mt*sf*rw)+cf*(-mf*cf*rw-mt*cf*rw))-rf*(mt*cf*rw+mf*cf*rw) -Cw-At-rf*(mf*rf+mt*r)-rt*mt*r-rw*cf*(mf*rf+mt*r)      0 ];
      
b = zeros(5,1);
b(1) = -V(t)+Ct*(-dphi*sf*dpsif*ct-dphi*cf*st*dtheta-cf*dpsif*dtheta);
b(2) = -U(t)+Cw*dphi*ct*dtheta-(-dphi*cf*ct+sf*dtheta)*Bf*(dphi*sf*ct+cf*dtheta)+(dphi*sf*ct+cf*dtheta)*Af*(-dphi*cf*ct+sf*dtheta)+At*dphi*ct*dtheta-(dphi*sf*ct+cf*dtheta)*Ct*(dphi*cf*ct-sf*dtheta+dpsit)+(dphi*cf*ct-sf*dtheta)*At*(dphi*sf*ct+cf*dtheta)-rf*(-mf*g*sf*ct-mf*(sf*dpsif*(-dphi*st+dpsiw)*rw+cf*dphi*ct*dtheta*rw-(-dphi*cf*ct+sf*dtheta)*(dtheta*rw+(dphi*sf*ct+cf*dtheta)*rf)+dphi*ct*dtheta*rf-(-dphi*st+dpsif)*sf*(-dphi*st+dpsiw)*rw)-mt*g*sf*ct-mt*(sf*dpsif*(-dphi*st+dpsiw)*rw+cf*dphi*ct*dtheta*rw-(-dphi*st+dpsif)*sf*(-dphi*st+dpsiw)*rw+dphi*ct*dtheta*(rf+rt)+(dphi*cf*ct-sf*dtheta)*(dtheta*rw+(dphi*sf*ct+cf*dtheta)*(rf+rt))))-rt*(-mt*g*sf*ct-mt*(sf*dpsif*(-dphi*st+dpsiw)*rw+cf*dphi*ct*dtheta*rw-(-dphi*st+dpsif)*sf*(-dphi*st+dpsiw)*rw+dphi*ct*dtheta*(rf+rt)+(dphi*cf*ct-sf*dtheta)*(dtheta*rw+(dphi*sf*ct+cf*dtheta)*(rf+rt))));
b(3) = -T*ct-2*dphi*st*Aw*dtheta-dtheta*Cw*(-dphi*st+dpsiw)+cf*(-Af*(dphi*sf*dpsif*ct+dphi*cf*st*dtheta+cf*dpsif*dtheta)-(dphi*sf*ct+cf*dtheta)*Cf*(-dphi*st+dpsif)+(-dphi*st+dpsif)*Bf*(dphi*sf*ct+cf*dtheta)+Ct*(-dphi*sf*dpsif*ct-dphi*cf*st*dtheta-cf*dpsif*dtheta))-sf*(-Bf*(dphi*cf*dpsif*ct-dphi*sf*st*dtheta-dpsif*sf*dtheta)-(-dphi*st+dpsif)*Af*(-dphi*cf*ct+sf*dtheta)+(-dphi*cf*ct+sf*dtheta)*Cf*(-dphi*st+dpsif)-At*(dphi*cf*dpsif*ct-dphi*sf*st*dtheta-dpsif*sf*dtheta)-(dphi*cf*ct-sf*dtheta)*At*(-dphi*st+dpsif)+(-dphi*st+dpsif)*Ct*(dphi*cf*ct-sf*dtheta+dpsit)+rf*(mf*g*st-mf*((dphi*sf*ct+cf*dtheta)*sf*(-dphi*st+dpsiw)*rw+(dphi*cf*dpsif*ct-dphi*sf*st*dtheta-dpsif*sf*dtheta)*rf+(-dphi*cf*ct+sf*dtheta)*(-cf*(-dphi*st+dpsiw)*rw-(-dphi*st+dpsif)*rf))+mt*g*st-mt*(-(dphi*cf*ct-sf*dtheta)*(-cf*(-dphi*st+dpsiw)*rw-(-dphi*st+dpsif)*(rf+rt))+(dphi*cf*dpsif*ct-dphi*sf*st*dtheta-dpsif*sf*dtheta)*(rf+rt)+(dphi*sf*ct+cf*dtheta)*sf*(-dphi*st+dpsiw)*rw))+rt*(mt*g*st-mt*(-(dphi*cf*ct-sf*dtheta)*(-cf*(-dphi*st+dpsiw)*rw-(-dphi*st+dpsif)*(rf+rt))+(dphi*cf*dpsif*ct-dphi*sf*st*dtheta-dpsif*sf*dtheta)*(rf+rt)+(dphi*sf*ct+cf*dtheta)*sf*(-dphi*st+dpsiw)*rw)));
b(4) = -dphi^2*st*Aw*ct-dphi*ct*Cw*(-dphi*st+dpsiw)-rw*(mw*dphi*ct*(-dphi*st+dpsiw)*rw-mt*g*st-mw*g*st+mf*((dphi*sf*ct+cf*dtheta)*sf*(-dphi*st+dpsiw)*rw+(dphi*cf*dpsif*ct-dphi*sf*st*dtheta-dpsif*sf*dtheta)*rf+(-dphi*cf*ct+sf*dtheta)*(-cf*(-dphi*st+dpsiw)*rw-(-dphi*st+dpsif)*rf))-mf*g*st+mt*(-(dphi*cf*ct-sf*dtheta)*(-cf*(-dphi*st+dpsiw)*rw-(-dphi*st+dpsif)*(rf+rt))+(dphi*cf*dpsif*ct-dphi*sf*st*dtheta-dpsif*sf*dtheta)*(rf+rt)+(dphi*sf*ct+cf*dtheta)*sf*(-dphi*st+dpsiw)*rw))+sf*(-Af*(dphi*sf*dpsif*ct+dphi*cf*st*dtheta+cf*dpsif*dtheta)-(dphi*sf*ct+cf*dtheta)*Cf*(-dphi*st+dpsif)+(-dphi*st+dpsif)*Bf*(dphi*sf*ct+cf*dtheta)+Ct*(-dphi*sf*dpsif*ct-dphi*cf*st*dtheta-cf*dpsif*dtheta))+cf*(-Bf*(dphi*cf*dpsif*ct-dphi*sf*st*dtheta-dpsif*sf*dtheta)-(-dphi*st+dpsif)*Af*(-dphi*cf*ct+sf*dtheta)+(-dphi*cf*ct+sf*dtheta)*Cf*(-dphi*st+dpsif)-At*(dphi*cf*dpsif*ct-dphi*sf*st*dtheta-dpsif*sf*dtheta)-(dphi*cf*ct-sf*dtheta)*At*(-dphi*st+dpsif)+(-dphi*st+dpsif)*Ct*(dphi*cf*ct-sf*dtheta+dpsit)+rf*(mf*g*st-mf*((dphi*sf*ct+cf*dtheta)*sf*(-dphi*st+dpsiw)*rw+(dphi*cf*dpsif*ct-dphi*sf*st*dtheta-dpsif*sf*dtheta)*rf+(-dphi*cf*ct+sf*dtheta)*(-cf*(-dphi*st+dpsiw)*rw-(-dphi*st+dpsif)*rf))+mt*g*st-mt*(-(dphi*cf*ct-sf*dtheta)*(-cf*(-dphi*st+dpsiw)*rw-(-dphi*st+dpsif)*(rf+rt))+(dphi*cf*dpsif*ct-dphi*sf*st*dtheta-dpsif*sf*dtheta)*(rf+rt)+(dphi*sf*ct+cf*dtheta)*sf*(-dphi*st+dpsiw)*rw))+rt*(mt*g*st-mt*(-(dphi*cf*ct-sf*dtheta)*(-cf*(-dphi*st+dpsiw)*rw-(-dphi*st+dpsif)*(rf+rt))+(dphi*cf*dpsif*ct-dphi*sf*st*dtheta-dpsif*sf*dtheta)*(rf+rt)+(dphi*sf*ct+cf*dtheta)*sf*(-dphi*st+dpsiw)*rw)));
b(5) = -T*st+2*Cw*dphi*ct*dtheta+(dphi*sf*ct+cf*dtheta)*Af*(-dphi*cf*ct+sf*dtheta)-rt*(-mt*g*sf*ct-mt*(sf*dpsif*(-dphi*st+dpsiw)*rw+cf*dphi*ct*dtheta*rw-(-dphi*st+dpsif)*sf*(-dphi*st+dpsiw)*rw+dphi*ct*dtheta*(rf+rt)+(dphi*cf*ct-sf*dtheta)*(dtheta*rw+(dphi*sf*ct+cf*dtheta)*(rf+rt))))-(dphi*sf*ct+cf*dtheta)*Ct*(dphi*cf*ct-sf*dtheta+dpsit)+At*dphi*ct*dtheta+rw*(2*mw*rw*dphi*ct*dtheta+sf*(mf*(-cf*dpsif*(-dphi*st+dpsiw)*rw+sf*dphi*ct*dtheta*rw+(dphi*sf*ct+cf*dtheta)*(dtheta*rw+(dphi*sf*ct+cf*dtheta)*rf)-(-dphi*st+dpsif)*(-cf*(-dphi*st+dpsiw)*rw-(-dphi*st+dpsif)*rf))-mf*g*cf*ct-mt*g*cf*ct-mt*(cf*dpsif*(-dphi*st+dpsiw)*rw-sf*dphi*ct*dtheta*rw+(-dphi*st+dpsif)*(-cf*(-dphi*st+dpsiw)*rw-(-dphi*st+dpsif)*(rf+rt))-(dphi*sf*ct+cf*dtheta)*(dtheta*rw+(dphi*sf*ct+cf*dtheta)*(rf+rt))))+cf*(mf*(sf*dpsif*(-dphi*st+dpsiw)*rw+cf*dphi*ct*dtheta*rw-(-dphi*cf*ct+sf*dtheta)*(dtheta*rw+(dphi*sf*ct+cf*dtheta)*rf)+dphi*ct*dtheta*rf-(-dphi*st+dpsif)*sf*(-dphi*st+dpsiw)*rw)+mt*g*sf*ct+mf*g*sf*ct+mt*(sf*dpsif*(-dphi*st+dpsiw)*rw+cf*dphi*ct*dtheta*rw-(-dphi*st+dpsif)*sf*(-dphi*st+dpsiw)*rw+dphi*ct*dtheta*(rf+rt)+(dphi*cf*ct-sf*dtheta)*(dtheta*rw+(dphi*sf*ct+cf*dtheta)*(rf+rt)))))+(dphi*cf*ct-sf*dtheta)*At*(dphi*sf*ct+cf*dtheta)-rf*(-mt*g*sf*ct-mt*(sf*dpsif*(-dphi*st+dpsiw)*rw+cf*dphi*ct*dtheta*rw-(-dphi*st+dpsif)*sf*(-dphi*st+dpsiw)*rw+dphi*ct*dtheta*(rf+rt)+(dphi*cf*ct-sf*dtheta)*(dtheta*rw+(dphi*sf*ct+cf*dtheta)*(rf+rt)))-mf*g*sf*ct-mf*(sf*dpsif*(-dphi*st+dpsiw)*rw+cf*dphi*ct*dtheta*rw-(-dphi*cf*ct+sf*dtheta)*(dtheta*rw+(dphi*sf*ct+cf*dtheta)*rf)+dphi*ct*dtheta*rf-(-dphi*st+dpsif)*sf*(-dphi*st+dpsiw)*rw))-(-dphi*cf*ct+sf*dtheta)*Bf*(dphi*sf*ct+cf*dtheta);

dz = zeros(12,1);
dz(1:5) = -A\b;
dz(6) = rw*cos(phi)*dpsiw;
dz(7) = rw*sin(phi)*dpsiw;
dz(8:12) = z(1:5);
