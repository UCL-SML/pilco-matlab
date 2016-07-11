% Augment a Gaussian with e*sin(x(i)), where i is a (possibly
% empty) set of I indices. The optional e scaling factor is a vector of length
% I. Optionally, compute derivatives of the parameters of the new Gaussian.
%
% Copyright (C) 2007, 2008 & 2009 by Carl Edward Rasmussen, 2009-07-07.

function [m, v, dmda, dmdb, dvda, dvdb] = trigSquash(a, b, i, e)

% a     is the mean vector of length "d" of a Gaussian distribution
% b     is the "d" by "d" covariance matrix
% i     is a vector of length I of indices of elements to augment
% e     is an optional scale vector of length I (defaults to unity)
%
% m     is a vector of length "D" of means
% v     is a "D" by "D" covariance matrix
% dmda  is a "D" by "d" matrix of derivatives of "m" wrt "a"
% dmdb  is a "D" by "d" by "d" matrix of derivatives of "m" wrt "b"
% dvda  is a "D" by "D" by "d" matrix of derivatives of "v" wrt "a"
% dvdb  is a "D" by "D" by "d" by "d" matrix of derivatives of "v" wrt "b"
%
% where the size "D" is "d+I".

d = length(a); I = length(i); D = d+I;
if nargin == 3, e = ones(I,1); else e = e(:); end          % unit column default
ai(1:I,1) = a(i); bi = b(i,i); bii(1:I,1) = diag(bi);      % short-hand notation

m = zeros(D,1); m(1:d) = a;                                    % first the means
m(d+1:end) = e.*exp(-bii/2).*sin(ai);

v = zeros(D); v(1:d,1:d) = b;                               % the covariance
v(d+1:end,1:d) = b(i,:).*repmat(e.*exp(-bii/2).*cos(ai),1,d);
v(1:d,d+1:end) = v(d+1:end,1:d)';                                % symmetric entries
for j=1:I
  for k=1:I
    if j==k
      q = e(j)*e(k)*(1-exp(-bii(j)))/2;
      v(d+j,d+j) = q*(1+exp(-bii(j))*cos(2*ai(j)));
    else
      q = e(j)*e(k)*exp(-(bii(j)+bii(k))/2)/2;

      % for numerical reasons:
      logq  = log(e(j)) + log(e(k)) -(bii(j)+bii(k))/2 - log(2);
      v(d+j,d+k) = (exp(logq + bi(j,k))-q)*cos(ai(j)-ai(k))...
        - (exp(logq-bi(j,k))-q)*cos(ai(j)+ai(k));
    end
  end
end

if nargout > 2                                            % compute derivatives?

  dmda = zeros(D,d); dmda(1:d,:) = eye(d);
  dmda(d+1:end,i) = diag(e.*exp(-bii/2).*cos(ai));

  dmdb = zeros(D,d,d);
  for j=1:I
    dmdb(d+j,i(j),i(j)) = -m(d+j)/2;
  end

  dvda = zeros(D,D,d);
  for j=1:I
    z = e(j)*b(:,i(j))*exp(-bii(j)/2);
    dvda(d+j,1:d,i(j)) = -z'*sin(ai(j));
    dvda(1:d,d+j,i(j)) = -z*sin(ai(j));
    for k=1:I
      if j==k
        z = e(j)*e(k)*(1-exp(-bii(j)))*exp(-bii(j));
        dvda(d+j,d+j,i(j)) = -z*sin(2*ai(j));
        
      else
        q = e(j)*e(k)*exp(-(bii(j)+bii(k))/2)/2;
        logq  = log(e(j)) + log(e(k)) -(bii(j)+bii(k))/2 - log(2);
        dvda(d+j,d+k,i(j)) = -(exp(logq + bi(j,k))-q)*sin(ai(j)-ai(k)) + ...
                                     (exp(logq - bi(j,k))-q)*sin(ai(j)+ai(k));
        dvda(d+j,d+k,i(k)) = (exp(logq + bi(j,k))-q)*sin(ai(j)-ai(k)) + ...
                                     (exp(logq - bi(j,k))-q)*sin(ai(j)+ai(k));
      end
    end
  end

  dvdb = zeros(D,D,d,d); 
  for j=1:d, for k=1:d, dvdb(j,k,j,k) = 0.5; dvdb(j,k,k,j) = dvdb(j,k,k,j) + 0.5; end, end
  for j=1:I
    dvdb(d+j,1:d,i(j),i(j)) = -v(d+j,1:d)/2;
    dvdb(1:d,d+j,i(j),i(j)) = -v(d+j,1:d)'/2;

    for k=1:I
        % derivative of trig variance w.r.t. input variance
      if j==k
          % trig variance
        q = e(j)*e(k)*exp(-bii(j))/2;
        dvdb(d+j,d+j,i(j),i(j)) = q*(1+cos(2*ai(j))*(2*exp(-bii(j))-1));
                
      else
         % trig-trig covariance terms
        logq  = log(e(j)) + log(e(k)) -(bii(j)+bii(k))/2 - log(2);
        dvdb(d+j,d+k,i(j),i(k)) = 0.5*((exp(logq+bi(j,k))*cos(ai(j)-ai(k)) + ...
          exp(logq-bi(j,k))*cos(ai(j)+ai(k))));
        dvdb(d+k,d+j,i(j),i(k)) = dvdb(d+j,d+k,i(j),i(k));
        dvdb(d+j,d+k,i(j),i(j)) = -v(d+j,d+k)/2;
        dvdb(d+j,d+k,i(k),i(k)) = -v(d+j,d+k)/2;
      end
    end
    z = e(j)*exp(-bii(j)/2)/2; zz = e(j)*(1-bii(j)/2)*exp(-bii(j)/2);
    for k = 1:d
        % derivative of covariance of trig-nontrig w.r.t input variance
      if k == i(j)
        dvdb(k,d+j,k,k) = zz*cos(ai(j));
        dvdb(d+j,k,k,k) = zz*cos(ai(j));
      else
        dvdb(k,d+j,k,i(j)) = z*cos(ai(j));
        dvdb(d+j,k,k,i(j)) = z*cos(ai(j));
        dvdb(k,d+j,i(j),k) = z*cos(ai(j));
        dvdb(d+j,k,i(j),k) = z*cos(ai(j));
      end
    end
  end
end