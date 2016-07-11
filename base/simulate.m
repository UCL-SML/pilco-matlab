%% simulate.m
% *Summary:* Simulate dynamics using a given control scheme.
%
%    function next = simulate(x0, f, plant)
%
% *Input arguments:*
%
%  x0      start state (with additional control states if required)
%  f       the control setpoint for this time step
%  plant   plant structure
%    .dt        time discretization
%    .dynamics  system function
%    .ctrl      function defining control implementation
%                  @zoh - zero-order-hold control (ZOH)
%                  @foh - first-order-hold control (FOH)
%                         with optional rise time 0 < plant.tau <= dt
%                  @lag - lagged control with time constant 0 < plant.tau
%    .delay     continuous-time delay, in range [0 dt)
%
% *Output arguments:*
%
%  next    successor state (with additional control states if required)
%
%
% Copyright (C) 2008-2013 by 
% Marc Deisenroth, Andrew McHutchon, Joe Hall, and Carl Edward Rasmussen.
%
% Last modification: 2012-06-30
%
%% High-Level Steps
% For each time step
% # Set up the control function
% # Simulate the dynamics (by calling ODE45)
% # Update control part of the state

function next = simulate(x0, f, plant)
%% Code


OPTIONS = odeset('RelTol', 1e-12, 'AbsTol', 1e-12);    % accuracy of ode45

x0 = x0(:); f = f(:); nU = length(f);
dt = plant.dt; dynamics = plant.dynamics;
if isfield(plant,'delay'), delay = plant.delay; else delay = 0; end
if isfield(plant,'tau'), tau = plant.tau; else tau = dt; end
par.dt = dt; par.delay = delay; par.tau = tau;

% 1. Set up control function ------------------------------------------------
% f{t} = control setpoint over time t to d+dt (determined by policy)
% u{t} = control currently being applied at time t
con = functions(plant.ctrl); con = con.function;
if (strcmp(con,'zoh') && delay==0)                            % U = [f{t}]
  x0s = x0;U = f; id = 0;
elseif strcmp(con,'zoh') || ...                          % U = [u{t} f{t}]
   (strcmp(con,'foh') && tau+delay<=dt) || (strcmp(con,'lag') && delay==0)
  x0s = x0(1:end-nU); U = [x0(end-nU+1:end) f]; id = 1;
else                                              % U = [f{t-1} u{t} f{t}]
  x0s=x0(1:end-2*nU); U=[reshape( x0(end-2*nU+1:end), [nU 2]) f]; id = 2;
end
ctrlfcn = str2func(con); u0 = cell(1,nU);        % create control function
for j = 1:nU, u0{j} = @(t)ctrlfcn(U(j,:),t,par); end

% 2. Simulate dynamics ------------------------------------------------------
[T y] = ode45(dynamics, [0 dt/2 dt], x0s, OPTIONS, u0{:});
x1 = y(3,:)';                                                 % next state 

% 3. Update control part of the state ---------------------------------------
udt = zeros(nU,1); for j=1:nU, udt(j) = u0{j}(dt); end
if id==0,     next =  x1;                         % return augmented state
elseif id==1, next = [x1; udt];
else          next = [x1; f; udt];
end


function u = zoh(f, t, par) % **************************** zero-order hold
d = par.delay;
if d==0
                  u = f;
else
  e = d/100; t0=t-(d-e/2);
  if t<d-e/2,     u=f(1);
  elseif t<d+e/2, u=(1-t0/e)*f(1) + t0/e*f(2);    % prevents ODE stiffness
  else            u=f(2);
  end
end

function u = foh(f, t, par) % *************************** first-order hold
d = par.delay; tau = par.tau; dt = par.dt;
if tau + d < dt
  t0=t-d;
  if t<d,         u=f(1);
  elseif t<tau+d, u=(1-t0/tau)*f(1) + t0/tau*f(2);
  else            u=f(2);
  end
else
  bit = d-(dt-tau);
  if t<bit,       u=(1-t/bit)*f(2) + t/tau*f(1);
  elseif t<d,     u=f(1);
  else t0=t+d;    u=(1-t0/tau)*f(1) + t0/tau*f(3);
  end
end

function u = lag(f, t, par) % **************************** first-order lag
d = par.delay; tau = par.tau;
if d==0
                  u = f(1) + (f(2)-f(1))*exp(-t/tau);
else
  bit = f(2) + (f(1)-f(2))*exp(-d/tau);
  if t<d,         u = f(2) + (f(1)-f(2))*exp(-t/tau);
  else            u = bit  + (f(3)-bit )*exp(-t/tau);
  end
end