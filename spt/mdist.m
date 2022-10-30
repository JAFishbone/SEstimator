% Copyright (C) 2020-2022 Justin A Fishbone
function d = mdist(X,mu,sigma)
% function d = mdist(X,mu,sigma)
% Squared Malanobis distance = (x-mu)'*sigma^-1*(x-mu)
%

cent  = X-mu;
dComp = (sigma\cent).*conj(cent);
d     = sum(dComp,1);
d     = real(d); % Force result to be real just in case there's a residual imaginary component due to numerical issues
end