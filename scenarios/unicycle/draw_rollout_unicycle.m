%% draw_rollout_unicycle.m
% *Summary:* Script to draw a rollout of the unicycle
%
% Copyright (C) 2008-2013 by
% Marc Deisenroth, Andrew McHutchon, Joe Hall, and Carl Edward Rasmussen.
%
% Last modified: 2013-04-04
%
%% High-Level Steps
% # For each time step, plot the observed trajectory, the applied torques,
% and the incurred cost

%% Code

if j > J
  
  draw_unicycle(xx, plant, plant.dt/10, cost, ...
    ['trial # ' num2str(j+J) ', T=' num2str(H*dt) ' sec'], ...
    ['total experience (after this trial): ' num2str(dt*size(x,1)) ...
    ' sec']);  
else
  draw_unicycle(xx, plant, plant.dt/10, cost, ...
    ['(random) trial # ' num2str(jj) ', T=' num2str(H*dt) ' sec'], ...
    ['total experience (after this trial): ' num2str(dt*size(x,1)) ...
    ' sec']);
end
