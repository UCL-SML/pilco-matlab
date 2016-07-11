%% settings_unicycle.m
% *Summary:* Script to set up the unicycle scenario
%
% Copyright (C) 2008-2013 by 
% Marc Deisenroth, Andrew McHutchon, Joe Hall, and Carl Edward Rasmussen.
%
% Last modified: 2013-04-02
%
%% High-Level Steps
% # Define state and important indices
% # Set up scenario
% # Set up the plant structure
% # Set up the policy structure
% # Set up the cost structure
% # Set up the GP dynamics model structure
% # Parameters for policy optimization
% # Plotting verbosity
% # Some array initializations


%% Code

warning('off','all'); format short; format compact

% include some paths
try
  rd = '../../';
  addpath([rd 'base'],[rd 'util'],[rd 'gp'],[rd 'control'],[rd 'loss']);
catch
end

rand('seed',1); randn('seed',1);

% 1. Define state and important indices

% 1a. Full state representation (including all augmentations)
% - non-angle velocities
% - angular velocities
% - non-angles
% - angles
% - controls
%
%  1  dx      x velocity
%  2  dy      y velocity
%  3  dxc     x velocity of origin (unicycle coordinates)
%  4  dyc     y velocity of origin (unicycle coordinates)
%  5  dtheta  roll angular velocity
%  6  dphi    yaw angular velocity
%  7  dpsiw   wheel angular velocity
%  8  dpsif   pitch angular velocity
%  9  dpsit   turn table angular velocity
% 10  x       x position
% 11  y       y position
% 12  xc      x position of origin (unicycle coordinates)
% 13  yc      y position of origin (unicycle coordinates)
% 14  theta   roll angle
% 15  phi     yaw angle
% 16  psiw    wheel angle
% 17  psif    pitch angle
% 18  psit    turn table angle
% 19  ct      control torque for turn table
% 20  cw      control torque for wheel

% 1b. Important indices
% odei  indicies for the ode solver
% augi  indicies for variables augmented to the ode variables
% dyno  indicies for the output from the dynamics model and indicies to loss
% angi  indicies for variables treated as angles (using sin/cos representation)
% dyni  indicies for inputs to the dynamics model
% poli  indicies for variables that serve as inputs to the policy
% difi  indicies for training targets that are differences (rather than values)

odei = [5 6 7 8 9 10 11 14 15 16 17 18];
augi = [1 2 3 4 12 13];
dyno = [5 6 7 8 9 12 13 14 15 17];
angi = [];
dyni = [1 2 3 4 5 6 7 8 10];
poli = [1 2 3 4 5 6 7 8 9 10];
difi = [1 2 3 4 5 6 7 8 9 10]; 


% 2. Set up the scenario
dt = 0.15;                    % [s] sampling time
T = 10.0;                     % [s] initial prediction horizon time
H = ceil(T/dt);               % prediction steps (optimization horizon)
maxH = ceil(10.0/dt);         % max pred horizon
s = [0.02 0.02 0.02 0.02 0.02 0.1 0.1 0.02 0.02 0.02 0.02 0.02].^2;
S0 = diag(s);                 % initial state variance, 95% is +/- 11.4 degrees
mu0 = zeros(length(odei),1);  % initial state mean
N = 40;                       % number controller optimizations
J = 10;                       % initial J trajectories of length H
K = 1;                        % number of initial states for which we optimize

% 3. Set up the plant structure
plant.dynamics = @dynamics_unicycle; % dynamics ODE function
plant.augment = @augment_unicycle; % function to augment the state ODE variables
plant.constraint = inline('abs(x(8))>pi/2 | abs(x(11))>pi/2');  % ODE constraint
plant.noise = diag([0.01*ones(1,5) 0.003*ones(1,7)].^2);     % measurement noise
plant.dt = dt;
plant.ctrl = @zoh;                 % controller is zero order hold
plant.odei = odei;                 % indices to the varibles for the ODE solver
plant.augi = augi;                 % indices of augmented variables
plant.angi = angi;
plant.poli = poli;
plant.dyno = dyno;
plant.dyni = dyni;
plant.difi = difi;
plant.prop = @propagated;

% 4. Set up the policy structure
policy.fcn = @(policy,m,s)conCat(@conlin,@gSat,policy,m,s); % linear policy
policy.maxU = [10 50];                                      % max. amplitude of 
policy.p.w = 1e-2*randn(length(policy.maxU),length(poli));  % weight matrix
policy.p.b = zeros(length(policy.maxU),1);                  % bias
                                                            % torques

% 5. Set up the cost structure
cost.fcn = @loss_unicycle;                  % cost function
cost.gamma = 1;                             % discount factor
cost.p = [0.22 0.81];                       % radius of wheel and length of rod
cost.width = 1;                           % cost function width
cost.expl = 0;                              % exploration parameter (UCB)

% 6. Set up the GP dynamics model structure
dynmodel.fcn = @gp1d;                % function for GP predictions
dynmodel.train = @train;             % function to train dynamics model
dynmodel.induce = zeros(300,0,1);    % shared inducing inputs (sparse GP)



% 7. Parameters for policy optimization
opt.length = -150;                   % max. number of function evaluations
opt.MFEPLS = 20;                     % max. number of function evaluations
                                     % per line search
opt.verbosity = 2;                   % verbosity: specifies how much 
                                     % information is displayed during
                                     % policy learning. Options: 0-3
opt.method = 'BFGS';                 % optimization algorithm
trainOpt = [300 300];                % defines the max. number of line searches
                                     % when training the GP dynamics models
                                     % trainOpt(1): full GP,
                                     % trainOpt(2): sparse GP (FITC)
                                     

% 8. Plotting verbosity
plotting.verbosity = 2;              % 0: no plots
                                     % 1: some plots
                                     % 2: all plots

% 9. Some initializations
x = []; y = [];
fantasy.mean = cell(1,N); fantasy.std = cell(1,N);
realCost = cell(1,N); M = cell(N,1); Sigma = cell(N,1);