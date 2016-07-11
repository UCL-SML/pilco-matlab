%% lossHinge.m
% *Summary:* Function to compute the moments and derivatives of the loss of a 
% Gaussian distributed point under a double hinge loss function. The loss 
% function has slope -/+a and corners b1 and b2. The function also calculates 
% derivatives of the loss w.r.t. the state distribution.
%
% Graph:
%          \                   /
%           \                 /
%            \               /
%             \_____________/
%             b1           b2
%
% To use a single hinge b1 or b2 can be set to -Inf or +Inf respectively.
%
% Note, this function is only analytic for 1D inputs. To apply this loss
% function to multiple variables, use the lossAdd function.
%
%
%   function [L dLdm dLds S dSdm dSds C dCdm dCds dLdb] = lossHinge(cost, m, s)
%
%
% *Input arguments:*
%
%   cost
%     .fcn      @lossHinge - called to get here
%     .a        slope of loss function
%     .b        corner points of loss function                     [1   x    2 ]
%   m             input mean                                       [D   x    1 ]
%   S             input covariance matrix                          [D   x    D ]
%
% *Output arguments:*
%
%   L               expected loss                                  [1   x    1 ]
%   dLdm            derivative of L wrt input mean                 [1   x    D ]
%   dLds            derivative of L wrt input covariance           [1   x   D^2] 
%   S               variance of loss                               [1   x    1 ]
%   dSdm            derivative of S wrt input mean                 [1   x    D ]
%   dSds            derivative of S wrt input covariance           [1   x   D^2]    
%   C               inv(S) times input-output covariance           [D   x    1 ]   
%   dCdm            derivative of C wrt input mean                 [D   x    D ]  
%   dCds            derivative of C wrt input covariance           [D   x   D^2]  
%
% Copyright (C) 2008-2013 by
% Marc Deisenroth, Andrew McHutchon, Joe Hall, and Carl Edward Rasmussen.
%
% Last modified: 2013-03-06
%
%% High-Level Steps
% # Expected cost
% # Variance of cost
% # inv(s)* cov(x,L)

function [L dLdm dLds S dSdm dSds C dCdm dCds dLdb] = lossHinge(cost, m, s)
%% Code
D = length(m); 
if D > 1; 
    error(['lossHinge only defined for 1D inputs, use lossAdd to '...
                                     'concatenate multiple 1D loss functions']);
end

a = cost.a; 
b = cost.b(:)' - m(:)'; I = ~isinf(b); % centralize
eb = exp(-b.^2/2/s); erfb = erf(b/sqrt(2*s));
c = sqrt(s/pi/2);

% 1. Expected Loss
% int_{-inf}^{b1-m} -a*(x-b1+m)*N(0,S)  +  int_{b2-m}^inf a*(x-b2+m)*N(0,S)
L = a*(b/2.*erfb + c*eb + b.*[1,-1]/2);
L = sum(L(I));

if nargout > 1
    % Derivative w.r.t. m
    dLdb = a/2*(erfb + [1,-1]); 
    dLdm = sum(dLdb(I)*-1);

    % Derivative w.r.t. S
    dc = 1/(2*sqrt(2*pi*s));
    dLds = a*sum(eb)*dc;
end
    
% 2. Variance of Loss
if nargout > 3
    S = a^2*((b.^2+s).*(1+[1,-1].*erfb)/2 + [1,-1].*b*c.*eb);
    S = sum(S(I)) - L^2;

    erfbdm = -sqrt(2/pi/s)*eb; erfbds = -b.*eb/sqrt(2*pi*s^3);
    dSdm = a^2*(-b.*(1+[1,-1].*erfb) + (b.^2+s).*[1,-1].*erfbdm/2 + ...
                                                    [1,-1]*c.*eb.*(-1+b.^2/s)); 
    dSdm = sum(dSdm(I)) - 2*L*dLdm;
    dSds = a^2/2*((1+[1,-1].*erfb) + (b.^2+s).*[1,-1].*erfbds + ...
                                   [2,-2].*b*dc.*eb + [1,-1].*b.^3/s^2*c.*eb);
    dSds = sum(dSds(I)) - 2*L*dLds;
end

% 3. inv(s)* covariance between input and cost
if nargout > 6
   C = a*([-1 1] - erfb)/2;
   C = sum(C);
   
   dCdm = -a/2*sum(erfbdm);
   dCds = -a/2*sum(erfbds(I));
end

