%% draw_rollout_dp.m
% *Summary:* Script to draw a trajectory of theobserved double-pendulum system 
% and the predicted uncertainties around the tips of the pendulums
%
% Copyright (C) 2008-2013 by 
% Marc Deisenroth, Andrew McHutchon, Joe Hall, and Carl Edward Rasmussen.
%
% Last modified: 2013-03-27
%
%% High-Level Steps
% # For each time step, plot the observed trajectory and the predicted
% means and covariances of the Cartesian coordinates of the tips of both
% pendulums

%% Code
% Loop over states in trajectory
for r = 1:size(xx,1)
  cost.t = r;
  if exist('j','var') && ~isempty(M{j})
    draw_dp(latent{j}(r,3), latent{j}(r,4), latent{j}(r,end-1), ...
      latent{j}(r,end), cost,  ...
      ['trial # ' num2str(j+J) ', T=' num2str(H*dt) ' sec'], ...
      ['total experience (after this trial): ' num2str(dt*size(x,1)) ...
      ' sec'], M{j}(:,r), Sigma{j}(:,:,r));
  else
    draw_dp(latent{jj}(r,3), latent{jj}(r,4), latent{jj}(r,end-1), ...
      latent{jj}(r,end), cost,  ...
      ['(random) trial # ' num2str(1) ', T=' num2str(H*dt) ' sec'], ...
      ['total experience (after this trial): ' num2str(dt*size(x,1)) ...
      ' sec'])
  end
  pause(dt);
end
  