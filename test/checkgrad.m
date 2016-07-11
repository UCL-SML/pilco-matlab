%% checkgrad.m
% *Summary:* checkgrad checks the derivatives in a function, by comparing them 
% to finite differences approximations. The partial derivatives and the 
% approximation are printed and the norm of the difference divided by the 
% norm of the sum is returned as an indication of accuracy.
%
%    function [d dy dh] = checkgrad(f, X, e, varargin)
%
%
% *Input arguments:*
%
%		f          function handle to function that needs to be checked. 
%              The function f should be of the type [fX, dfX] = f(X, varargin)
%              where fX is the function value and dfX is the gradient of fX
%              with respect to the parameters X
%   X          parameters (can be a vector or a struct)
%   e          small perturbation used for finite differences (1e-4 is good)
%   varargin   other arguments that are passed on to the function f
%
%
% *Output arguments:*
%
%   d          relative error of analytical vs. finite difference gradient
%   dy         analytical gradient
%   dh         finite difference gradient
%
%
% Copyright (C) 2008-2013 by
% Marc Deisenroth, Andrew McHutchon, Joe Hall, and Carl Edward Rasmussen.
%
% Last modified: 2013-03-21
%
%% High-Level Steps
% # Analytical gradient
% # Numerical gradient via finite differences
% # Relative error

function [d dy dh] = checkgrad(f, X, e, varargin)
%% Code

% 1. Analytical gradient
Z = unwrap(X); NZ = length(Z);                 % number of input variables
[y dy] = feval(f, X, varargin{:});             % get the partial derivatives dy
[D E] = size(y); y = y(:); Ny = length(y);     % number of output variables
if iscell(dy) || isstruct(dy); dy = unwrap(dy); end;
dy = reshape(dy,Ny,NZ);

% 2. Finite difference approximation
dh = zeros(Ny,NZ);
for j = 1:NZ
  dx = zeros(length(Z),1);
  dx(j) = dx(j) + e;                               % perturb a single dimension
  y2 = feval(f, rewrap(X,Z+dx), varargin{:});                        
  y1 = feval(f, rewrap(X,Z-dx), varargin{:});
  dh(:,j) = (y2(:) - y1(:))/(2*e);
end

% 3. Compute error 
% norm of diff divided by norm of sum
d = sqrt(sum((dh-dy).^2,2)./sum((dh+dy).^2,2)); 
small = max(abs([dy dh]),[],2) < 1e-5; % small derivatives are poorly tested ...
d(d > 1e-3 & small) = NaN;             % ... by finite differences
d = reshape(d,D,E);

disp('   Analytic  Numerical');
for i=1:Ny;
    disp([dy(i,:)' dh(i,:)']);                           % print the two vectors
    fprintf('d = %e\n\n',d(i))
end 

if Ny > 1; disp('For all outputs, d = '); disp(d); end; fprintf('\n');
