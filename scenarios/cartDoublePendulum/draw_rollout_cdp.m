%% draw_rollout_cdp.m
% *Summary:* Script to draw a trajectory of the cart-double-pendulum system and
%  the predicted uncertainties around the tips of the pendulums
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
for r = 1:size(xx,1)
  if exist('j','var') && ~isempty(M{j})
    draw_cdp(latent{j}(r,1), latent{j}(r,5), latent{j}(r,6), ...
      latent{j}(r,end), cost, M{j}(:,r), Sigma{j}(:,:,r), ...
      ['trial # ' num2str(j+J) ', T=' num2str(H*dt) ' sec'], ...
      ['total experience (after this trial): ' num2str(dt*size(x,1)) ...
      ' sec']);
  else
    draw_cdp(latent{jj}(r,1), latent{jj}(r,5), latent{jj}(r,6), ...
      latent{jj}(r,end), cost, [], [], ...
      ['(random) trial # ' num2str(jj) ', T=' num2str(H*dt) ' sec'], ...
      ['total experience (after this trial): ' num2str(dt*size(x,1)) ...
      ' sec']);
  end
  pause(dt);
end
