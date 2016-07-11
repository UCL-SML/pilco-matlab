%% gp2d.m
% *Summary:* Compute joint predictions and derivatives for multiple GPs
% with uncertain inputs. Does not consider the uncertainty about the underlying
% function (in prediction), hence, only the GP mean function is considered.
% Therefore, this representation is equivalent to a regularized RBF
% network.
% If gpmodel.nigp exists, individial noise contributions are added.
%
%
%   function [M, S, V, dMdm, dSdm, dVdm, dMds, dSds, dVds, dMdi, dSdi, dVdi, ...
%                      dMdt, dSdt, dVdt, dMdX, dSdX, dVdX] = gp2d(gpmodel, m, s)
%
% *Input arguments:*
%
%   gpmodel    GP model struct
%     hyp      log-hyper-parameters                                  [D+2 x  E ]
%     inputs   training inputs                                       [ n  x  D ]
%     targets  training targets                                      [ n  x  E ]
%     nigp     (optional) individual noise variance terms            [ n  x  E ]
%   m          mean of the test distribution                         [ D  x  1 ]
%   s          covariance matrix of the test distribution            [ D  x  D ]
%
% *Output arguments:*
%
%   M          mean of pred. distribution                            [ E  x  1 ]
%   S          covariance of the pred. distribution                  [ E  x  E ]
%   V          inv(s) times covariance between input and output      [ D  x  E ]
%   dMdm       output mean by input mean                             [ E  x  D ]
%   dSdm       output covariance by input mean                       [E*E x  D ]
%   dVdm       inv(s)*input-output covariance by input mean          [D*E x  D ]
%   dMds       ouput mean by input covariance                        [ E  x D*D]
%   dSds       output covariance by input covariance                 [E*E x D*D]
%   dVds       inv(s)*input-output covariance by input covariance    [D*E x D*D]
%
%   dMdi      output mean by inputs                                  [ E  x n*D]
%   dSdi      output covariance by inputs                            [E*E x n*D]
%   dVdi      inv(s) times input-output covariance by inputs         [D*E x n*D]
%   dMdt      output mean by targets                                 [ E  x n*E]
%   dSdt      output covariance by targets                           [E*E x n*E]
%   dVdt      inv(s) times input-output covariance by targets        [D*E x n*E]
%   dMdX      output mean by hyperparameters                         [ E  x P*E]
%   dSdX      output covariance by hyperparameters                   [E*E x P*E]
%   dVdX      inv(s) times input-output covariance by hyper-par.     [D*E x P*E]
%
%
%
% Copyright (C) 2008-2013 by
% Marc Deisenroth, Andrew McHutchon, Joe Hall, and Carl Edward Rasmussen.
%
% Last modified: 2013-05-24
%
%% High-Level Steps
% # If necessary, re-compute cached variables
% # Compute predicted mean and inv(s) times input-output covariance
% # Compute predictive covariance matrix, non-central moments
% # Centralize moments
% # Vectorize derivatives

function [M, S, V, dMdm, dSdm, dVdm, dMds, dSds, dVds, dMdi, dSdi, dVdi, ...
  dMdt, dSdt, dVdt, dMdX, dSdX, dVdX] = gp2d(gpmodel, m, s)
%% Code
input = gpmodel.inputs;  target = gpmodel.targets; X = gpmodel.hyp;

if nargout < 4; [M, S, V] = gp2(gpmodel, m, s); return; end

persistent K iK oldX oldIn oldOut beta oldn;           % cache some variables
D = size(input,2);         % number of examples and dimension of input space
[n, E] = size(target);               % number of examples and number of outputs
X = reshape(X, D+2, E);

% 1) If necessary, re-compute cached variables
if numel(X) ~= numel(oldX) || isempty(iK) ||  n ~= oldn || ...
    sum(any(X ~= oldX)) || sum(any(oldIn ~= input)) || ...
    sum(any(oldOut ~= target))
  oldX = X; oldIn = input; oldOut = target; oldn = n;
  K = zeros(n,n,E); iK = K; beta = zeros(n,E);
  
  % compute K and inv(K) and beta
  for i=1:E
    inp = bsxfun(@rdivide,input,exp(X(1:D,i)'));
    K(:,:,i) = exp(2*X(D+1,i)-maha(inp,inp)/2);
    if isfield(gpmodel,'nigp')
      L = chol(K(:,:,i) + exp(2*X(D+2,i))*eye(n) + diag(gpmodel.nigp(:,i)))';
    else
      L = chol(K(:,:,i) + exp(2*X(D+2,i))*eye(n))';
    end
    iK(:,:,i) = L'\(L\eye(n));
    beta(:,i) = L'\(L\gpmodel.targets(:,i));
  end
  
end

% initializations
k = zeros(n,E); M = zeros(E,1); V = zeros(D,E); S = zeros(E);
dMds = zeros(E,D,D); dSdm = zeros(E,E,D); r = zeros(1,D);
dSds = zeros(E,E,D,D); dVds = zeros(D,E,D,D); T = zeros(D);
tlbdi = zeros(n,D); dMdi = zeros(E,n,D); dMdt = zeros(E,n,E);
dVdt = zeros(D,E,n,E); dVdi = zeros(D,E,n,D); dSdt = zeros(E,E,n,E);
dSdi = zeros(E,E,n,D); dMdX = zeros(E,D+2,E); dSdX = zeros(E,E,D+2,E);
dVdX = zeros(D,E,D+2,E); Z = zeros(n,D);
bdX = zeros(n,E,D); kdX = zeros(n,E,D+1);

% centralize training inputs
inp = bsxfun(@minus,input,m');

% 2) compute predicted mean and input-output covariance
for i = 1:E
  % first some useful intermediate terms
  K2 = K(:,:,i)+exp(2*X(D+2,i))*eye(n);                         % K + sigma^2*I
  inp2 = bsxfun(@rdivide,input,exp(X(1:D,i)'));
  ii = bsxfun(@rdivide,input,exp(2*X(1:D,i)'));
  R = s+diag(exp(2*X(1:D,i)));
  L = diag(exp(-X(1:D,i)));
  B = L*s*L+eye(D); iR = L/B*L;
  t = inp*iR;
  l = exp(-sum(t.*inp,2)/2); lb = l.*beta(:,i);
  tliK = t'*bsxfun(@times,l,iK(:,:,i));
  liK = K2\l;
  tlb = bsxfun(@times,t,lb);
  
  c = exp(2*X(D+1,i))/sqrt(det(R))*exp(sum(X(1:D,i)));
  detdX = diag(bsxfun(@times,det(R)*iR',2.*exp(2.*X(1:D,i))));    % d(det R)/dX
  cdX = -0.5*c/det(R).*detdX'+ c.*ones(1,D);       % derivs w.r.t length-scales
  dldX = bsxfun(@times,l,bsxfun(@times,t,2.*exp(2*X(1:D,i)')).*t./2);
  
  M(i) = sum(lb)*c;                                            % predicted mean
  
  iK2beta = K2\beta(:,i);
  dMds(i,:,:) = c*t'*tlb/2-iR*M(i)/2;
  dMdX(i,D+2,i) = -c*sum(l.*(2*exp(2*X(D+2,i))*iK2beta));           % OK
  dMdX(i,D+1,i) = -dMdX(i,(i-1)*(D+2)+D+2);
  
  dVdX(:,i,D+2,i) = -((l.*(2*exp(2*X(D+2,i))*iK2beta))'*t*c)';
  dVdX(:,i,D+1,i) = -dVdX(:,i,D+2,i);
  
  dsi = -bsxfun(@times,inp2,2.*inp2);                 % d(sum(inp2.*inp2,2))/dX
  dslb = zeros(1,D);
  
  for d = 1:D
    sqdi = K(:,:,i).*bsxfun(@minus,ii(:,d),ii(:,d)');
    sqdiBi = sqdi*beta(:,i);
    tlbdi(:,d) = sqdi*liK.*beta(:,i) + sqdiBi.*liK;
    tlbdi2 = -tliK*(-bsxfun(@times,sqdi,beta(:,i))'-diag(sqdiBi));
    dVdi(:,i,:,d) = c*(iR(:,d)*lb' - bsxfun(@times,t,tlb(:,d))' + tlbdi2);
    dsqdX = bsxfun(@plus,dsi(:,d),dsi(:,d)') + 4.*inp2(:,d)*inp2(:,d)';
    dKdX = -K(:,:,i).*dsqdX./2;                               % dK/dX(1:D)
    dKdXbeta = dKdX*beta(:,i);
    bdX(:,i,d) = -K2\dKdXbeta;                        % dbeta/dX
    dslb(d) = -liK'*dKdXbeta + beta(:,i)'*dldX(:,d);
    dlb = dldX(:,d).*beta(:,i) + l.*bdX(:,i,d);
    dtdX = inp*(-bsxfun(@times,iR(:,d),2.*exp(2*X(d,i))*iR(d,:)));
    dlbt = lb'*dtdX + dlb'*t;
    dVdX(:,i,d,i) = (dlbt'*c + cdX(d)*(lb'*t)');
  end % d
  
  dMdi(i,:,:) = c*(tlbdi - tlb);
  dMdt(i,:,i) = c*liK';
  dMdX(i,1:D,i) = cdX.*sum(beta(:,i).*l) + c.*dslb;
  v = bsxfun(@rdivide,inp,exp(X(1:D,i)'));
  k(:,i) = 2*X(D+1,i)-sum(v.*v,2)/2;
  V(:,i) = t'*lb*c;                                   % input-output covariance
  
  for d = 1:D
    dVds(d,i,:,:) = c*bsxfun(@times,t,t(:,d))'*tlb/2 - iR*V(d,i)/2 ...
      - V(:,i)*iR(d,:)/2 -iR(:,d)*V(:,i)'/2;
    kdX(:,i,d) = bsxfun(@times,v(:,d),v(:,d));
  end % d
  
  dVdt(:,i,:,i) = c*tliK;
  kdX(:,i,D+1) = 2*ones(1,n);                       % pre-computation for later
  
end % i
dMdm = V';                                                % derivatives w.r.t m
dVdm = 2*permute(dMds,[2 1 3]);


% 3) predictive covariance matrix (non-central moments)
for i = 1:E
  K2 = K(:,:,i)+exp(2*X(D+2,i))*eye(n);
  ii = bsxfun(@rdivide,inp,exp(2*X(1:D,i)'));
  
  for j = 1:i % if i==j: diagonal elements of S; see Marc's thesis around eq. (2.26)
    R = s*diag(exp(-2*X(1:D,i))+exp(-2*X(1:D,j)))+eye(D); t = 1./sqrt(det(R));
    if rcond(R) < 1e-15; fprintf('R-matrix in gp2d ill-conditioned'); keyboard; end
    iR = R\eye(D);
    ij = bsxfun(@rdivide,inp,exp(2*X(1:D,j)'));
    L = exp(bsxfun(@plus,k(:,i),k(:,j)')+maha(ii,-ij,R\s/2)); % called Q in thesis
    A = beta(:,i)*beta(:,j)'; A = A.*L; ssA = sum(sum(A));
    S(i,j) = t*ssA; S(j,i) = S(i,j);
    
    zzi = ii*(R\s);
    zzj = ij*(R\s);
    zi = ii/R; zj = ij/R;
    
    tdX  = -0.5*t*sum(iR'.*bsxfun(@times,s,-2*exp(-2*X(1:D,i)')-2*exp(-2*X(1:D,i)')));
    tdXi = -0.5*t*sum(iR'.*bsxfun(@times,s,-2*exp(-2*X(1:D,i)')));
    tdXj = -0.5*t*sum(iR'.*bsxfun(@times,s,-2*exp(-2*X(1:D,j)')));
    bLiKi = iK(:,:,j)*(L'*beta(:,i)); bLiKj = iK(:,:,i)*(L*beta(:,j));
    
    Q2 = R\s/2;
    aQ = ii*Q2; bQ = ij*Q2;
    
    for d = 1:D
      
      Z(:,d) = exp(-2*X(d,i))*(A*zzj(:,d) + sum(A,2).*(zzi(:,d) - inp(:,d)))...
        + exp(-2*X(d,j))*((zzi(:,d))'*A + sum(A,1).*(zzj(:,d) - inp(:,d))')';
      Q = bsxfun(@minus,inp(:,d),inp(:,d)');
      B = K(:,:,i).*Q;
      Z(:,d) = Z(:,d)+exp(-2*X(d,i))*(B*beta(:,i).*bLiKj+beta(:,i).*(B*bLiKj));
      
      if i~=j; B = K(:,:,j).*Q; end
      
      Z(:,d) = Z(:,d)+exp(-2*X(d,j))*(bLiKi.*(B*beta(:,j))+B*bLiKi.*beta(:,j));
      B = bsxfun(@plus,zi(:,d),zj(:,d)').*A;
      r(d) = sum(sum(B))*t;
      T(d,1:d) = sum(zi(:,1:d)'*B,2) + sum(B*zj(:,1:d))';
      T(1:d,d) = T(d,1:d)';
      
      if i==j
        RTi =  bsxfun(@times,s,(-2*exp(-2*X(1:D,i)')-2*exp(-2*X(1:D,j)')));
        diRi = -R\bsxfun(@times,RTi(:,d),iR(d,:));
      else
        RTi = bsxfun(@times,s,-2*exp(-2*X(1:D,i)'));
        RTj = bsxfun(@times,s,-2*exp(-2*X(1:D,j)'));
        diRi = -R\bsxfun(@times,RTi(:,d),iR(d,:));
        diRj = -R\bsxfun(@times,RTj(:,d),iR(d,:));
        QdXj = diRj*s/2; % dQ2/dXj
      end
      
      QdXi = diRi*s/2; % dQ2/dXj
      
      if i==j
        daQi = ii*QdXi + bsxfun(@times,-2*ii(:,d),Q2(d,:)); % d(ii*Q)/dXi
        dsaQi = sum(daQi.*ii,2) - 2.*aQ(:,d).*ii(:,d); dsaQj = dsaQi;
        dsbQi = dsaQi; dsbQj = dsbQi;
        dm2i = -2*daQi*ii' + 2*(bsxfun(@times,aQ(:,d),ii(:,d)')...
          +bsxfun(@times,ii(:,d),aQ(:,d)')); dm2j = dm2i; % -2*aQ*ij'/di
      else
        dbQi = ij*QdXi;  % d(ij*Q)/dXi
        dbQj = ij*QdXj + bsxfun(@times,-2*ij(:,d),Q2(d,:)); % d(ij*Q)/dXj
        daQi = ii*QdXi + bsxfun(@times,-2*ii(:,d),Q2(d,:)); % d(ii*Q)/dXi
        daQj = ii*QdXj; % d(ii*Q)/dXj
        
        dsaQi = sum(daQi.*ii,2) - 2.*aQ(:,d).*ii(:,d);
        dsaQj = sum(daQj.*ii,2);
        dsbQi = sum(dbQi.*ij,2);
        dsbQj = sum(dbQj.*ij,2) - 2.*bQ(:,d).*ij(:,d);
        dm2i = -2*daQi*ij'; % second part of the maha(..) function wrt Xi
        dm2j = -2*ii*(dbQj)'; % second part of the maha(..) function wrt Xj
      end
      
      dm1i = bsxfun(@plus,dsaQi,dsbQi'); % first part of the maha(..) function wrt Xi
      dm1j = bsxfun(@plus,dsaQj,dsbQj'); % first part of the maha(..) function wrt Xj
      dmahai = dm1i-dm2i;
      dmahaj = dm1j-dm2j;
      
      if i==j
        LdXi = L.*(dmahai + bsxfun(@plus,kdX(:,i,d),kdX(:,j,d)'));
        dSdX(i,i,d,i) = beta(:,i)'*LdXi*beta(:,j);
      else
        LdXi = L.*(dmahai + bsxfun(@plus,kdX(:,i,d),zeros(n,1)'));
        LdXj = L.*(dmahaj + bsxfun(@plus,zeros(n,1),kdX(:,j,d)'));
        dSdX(i,j,d,i) = beta(:,i)'*LdXi*beta(:,j);
        dSdX(i,j,d,j) = beta(:,i)'*LdXj*beta(:,j);
      end
      
    end % d
    
    if i==j
      dSdX(i,i,1:D,i) = reshape(dSdX(i,i,1:D,i),D,1) + reshape(bdX(:,i,:),n,D)'*(L+L')*beta(:,i);
      dSdX(i,i,1:D,i) = reshape(t*dSdX(i,i,1:D,i),D,1)' + tdX*ssA;
      dSdX(i,i,D+2,i) = 2*exp(2*X(D+2,i))*t*(-sum(beta(:,i).*bLiKi)-sum(beta(:,i).*bLiKi));
    else
      dSdX(i,j,1:D,i) = reshape(dSdX(i,j,1:D,i),D,1) + reshape(bdX(:,i,:),n,D)'*(L*beta(:,j));
      dSdX(i,j,1:D,j) = reshape(dSdX(i,j,1:D,j),D,1) + reshape(bdX(:,j,:),n,D)'*(L'*beta(:,i));
      dSdX(i,j,1:D,i) = reshape(t*dSdX(i,j,1:D,i),D,1)' + tdXi*ssA;
      dSdX(i,j,1:D,j) = reshape(t*dSdX(i,j,1:D,j),D,1)' + tdXj*ssA;
      dSdX(i,j,D+2,i) = 2*exp(2*X(D+2,i))*t*(-beta(:,i)'*bLiKj);
      dSdX(i,j,D+2,j) = 2*exp(2*X(D+2,j))*t*(-beta(:,j)'*bLiKi);
    end
    
    dSdm(i,j,:) = r - M(i)*dMdm(j,:)-M(j)*dMdm(i,:); dSdm(j,i,:) = dSdm(i,j,:);
    T = (t*T-S(i,j)*diag(exp(-2*X(1:D,i))+exp(-2*X(1:D,j)))/R)/2;
    T = T - reshape(M(i)*dMds(j,:,:) + M(j)*dMds(i,:,:),D,D);
    dSds(i,j,:,:) = T; dSds(j,i,:,:) = T;
    
    if i==j
      dSdt(i,i,:,i) = (beta(:,i)'*(L+L'))/(K2)*t ...
        - 2*dMdt(i,:,i)*M(i);
      dSdX(i,j,:,i) = reshape(dSdX(i,j,:,i),1,D+2) - M(i)*dMdX(j,:,j)-M(j)*dMdX(i,:,i);
    else
      dSdt(i,j,:,i) = (beta(:,j)'*L')/(K2)*t ...
        - dMdt(i,:,i)*M(j);
      dSdt(i,j,:,j) = beta(:,i)'*L/(K(:,:,j)+exp(2*X(D+2,j))*eye(n))*t ...
        - dMdt(j,:,j)*M(i);
      dSdt(j,i,:,:) = dSdt(i,j,:,:);
      dSdX(i,j,:,j) = reshape(dSdX(i,j,:,j),1,D+2) - M(i)*dMdX(j,:,j);
      dSdX(i,j,:,i) = reshape(dSdX(i,j,:,i),1,D+2) - M(j)*dMdX(i,:,i);
    end
    
    dSdi(i,j,:,:) = Z*t - reshape(M(i)*dMdi(j,:,:) + dMdi(i,:,:)*M(j),n,D);
    dSdi(j,i,:,:) = dSdi(i,j,:,:);
    dSdX(j,i,:,:) = dSdX(i,j,:,:);
  end % j
  
  S(i,i) = S(i,i) + 1e-06;    % add small diagonal jitter for numerical reasons
end % i

dSdX(:,:,D+1,:) = -dSdX(:,:,D+2,:);
dSdX(:,:,D+1,:) = -dSdX(:,:,D+2,:);

% 4) centralize moments
S = S - M*M';
%S(diag(S)<0,diag(S)<0) = 1e-6;

% 5) Vectorize derivatives
dMds=reshape(dMds,[E D*D]);
dSdm=reshape(dSdm,[E*E D]); dSds=reshape(dSds,[E*E D*D]);
dVdm=reshape(dVdm,[D*E D]); dVds=reshape(dVds,[D*E D*D]);
dMdi=reshape(dMdi,E,[]);  dMdt=reshape(dMdt,E,[]);  dMdX=reshape(dMdX,E,[]);
dSdi=reshape(dSdi,E*E,[]);dSdt=reshape(dSdt,E*E,[]);dSdX=reshape(dSdX,E*E,[]);
dVdi=reshape(dVdi,D*E,[]);dVdt=reshape(dVdt,D*E,[]);dVdX=reshape(dVdX,D*E,[]);