%% fitc.m
% *Summary:* Compute the FITC negative log marginal likelihood and its
% derivatives with  respect to the inducing inputs (we don't compute the
% derivatives with respect to the GP hyper-parameters)
%
%   function [nml dnml] = fitc(induce, gpmodel)
%
% *Input arguments:*
%
%   induce          matrix of inducing inputs                       [M x D x uE]
%                   M: number of inducing inputs
%                   E: either 1 (inducing inputs are shared across target dim.)
%                      or     E (different inducing inputs for each target dim.)
%   gpmodel         GP structure
%     .hyp          log-hyper-parameters                               [D+2 x E]
%     .inputs       training inputs                                    [N   x D]
%     .targets      training targets                                   [N   x E]
%     .noise (opt)  noise
%
%
% *Output arguments:*
%
%   nlml             negative log-marginal likelihood
%   dnlml            derivative of negative log-marginal likelihood wrt
%                    inducing inputs
%
% Adapted from Ed Snelson's SPGP code.
%
% Copyright (C) 2008-2013 by
% Marc Deisenroth, Andrew McHutchon, Joe Hall, and Carl Edward Rasmussen.
%
% Last modified: 2013-05-21
%
%% High-Level Steps
%
% # Compute FITC marginal likelihood 
% # Compute corresponding gradients wrt the pseudo inputs

function [nlml dnlml] = fitc(induce, gpmodel)
%% Code

ridge = 1e-06;                       % jitter to make matrix better conditioned

[N D] = size(gpmodel.inputs); E = size(gpmodel.targets,2);
[M uD uE] = size(induce);
if uD ~= D || (uE~=1 && uE ~= E); error('Wrong size of inducing inputs'); end

nlml = 0; dfxb = zeros(M, D); dnlml = zeros(M, D, E); % zero and allocate outputs

for j = 1:E
  if uE > 1; u = induce(:,:,j); else u = induce; end
  b = exp(gpmodel.hyp(1:D,j));                                 % length-scales
  c = gpmodel.hyp(D+1,j);                                 % log signal std dev
  sig = exp(2.*gpmodel.hyp(D+2,j));                           % noise variance
  
  xb = bsxfun(@rdivide,u,b');                 % divide inducing by lengthscales
  x = bsxfun(@rdivide,gpmodel.inputs,b');     % divide inputs by length-scales
  y = gpmodel.targets(:,j);                                  % training targets
  
  Kmm = exp(2*c-maha(xb,xb)/2) + ridge*eye(M);
  Kmn = exp(2*c-maha(xb,x)/2);
  
  % Check whether Kmm is no longer positive definite. If so, return
  try 
    L = chol(Kmm)'; 
  catch
    nlml = Inf; dnlml = zeros(size(params)); 
    return;
  end
  V = L\Kmn;                                               % inv(sqrt(Kmm))*Kmn
  
  if isfield(gpmodel,'noise')
    Gamma = 1 + (exp(2*c)-sum(V.^2)'+gpmodel.noise(:,j))/sig;
  else
    Gamma = 1 + (exp(2*c)-sum(V.^2)')/sig;      % Gamma = diag(Knn-Qnn)/sig + I
  end
  
  V = bsxfun(@rdivide,V,sqrt(Gamma)');  % inv(sqrt(Kmm))*Kmn * inv(sqrt(Gamma))
  y = y./sqrt(Gamma);
  Am = chol(sig*eye(M) + V*V')';        % chol(inv(sqrt(Kmm))*A*inv(sqrt(Kmm)))
                % V*V' = inv(chol(Kmm)')*K*inv(diag(Gamma))*K'*inv(chol(Kmm)')'
  Vy = V*y;
  beta = Am\Vy;
  
  nlml = nlml + sum(log(diag(Am))) + (N-M)/2*log(sig) + sum(log(Gamma))/2 ...
         + (y'*y - beta'*beta)/2/sig + 0.5*N*log(2*pi);
  
  if nargout == 2               % ... and if requested, its partial derivatives
    
    At = L*Am; iAt = At\eye(M);              % chol(sig*B) [Ed's thesis, p. 40]
    iA = iAt'*iAt;                                                 % inv(sig*B)
    
    iAmV = Am\V;                                                    % inv(Am)*V
    B1 = At'\(iAmV);
    b1 = At'\beta;                                                  % b1 = B1*y
    
    iLV = L'\V;                                 % inv(Kmm)*Kmn*inv(sqrt(Gamma))
    iL = L\eye(M);
    iKmm = iL'*iL;
    
    mu = ((Am'\beta)'*V)';
    bs = y.*(beta'*iAmV)'/sig - sum(iAmV.*iAmV)'/2 - (y.^2+mu.^2)/2/sig + 0.5;
    TT = iLV*(bsxfun(@times,iLV',bs));
    Kmn = bsxfun(@rdivide,Kmn,sqrt(Gamma)');                    % overwrite Kmn
    
    for i = 1:D                               % derivatives wrt inducing inputs
      dsq_mm = bsxfun(@minus,xb(:,i),xb(:,i)').*Kmm;
      dsq_mn = bsxfun(@minus,-xb(:,i),-x(:,i)').*Kmn;
      dGamma = -2/sig*dsq_mn.*iLV;
      
      dfxb(:,i) = -b1.*(dsq_mn*(y-mu)/sig + dsq_mm*b1) + dGamma*bs ...
                  + sum((iKmm - iA*sig).*dsq_mm,2) - 2/sig*sum(dsq_mm.*TT,2);
      dsq_mn = dsq_mn.*B1;                                   % overwrite dsq_mn
      dfxb(:,i) = dfxb(:,i) + sum(dsq_mn,2);
      dfxb(:,i) = dfxb(:,i)/b(i);
    end
    
    dnlml(:,:,j) = dfxb;
  end
end
if 1 == uE; dnlml = sum(dnlml,3); end % combine derivatives if sharing inducing
