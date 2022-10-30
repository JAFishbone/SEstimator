% Copyright (C) 2020-2022 Justin A Fishbone
function [muHat, sigmaHat, omegaHat] = solveMLE(X, wFun, mu0, sigma0, tol)
% [muHat, sigmaHat, omegaHat] = solveMLE(X, wFun)
% [muHat, sigmaHat, omegaHat] = solveMLE(X, wFun, mu0, sigma0)
% [muHat, sigmaHat, omegaHat] = solveMLE(X, wFun, mu0, sigma0, tol)
% Solve MLE for location and scatter/shape for elliptically distributed r.v.
%
% INPUTS:
% X      - pxn random data
% wFun   - Function handle for MLE weight function -  -2*phi_prime(d)./phi(d)
% mu0    - (Optional) px1 initial location guess
% sigma0 - (Optional) pxp initial scatter guess
% tol    - (Optional) Convergence criteria
%
% OUTPUTS:
% muHat    - px1 location matrix estimate
% sigmaHat - pxp scatter matrix estimate
% omegaHat - pxp shape matrix estimate
% 

%--------------------------------------------------------------
%%% PROCESS INPUTS
%--------------------------------------------------------------

% Data size
p = size(X,1);
n = size(X,2);

% Numerical convergence criteria
if nargin < 5 || isempty(tol)
    tol = 1e-6;
end

% Initial estimates
if nargin < 3 || isempty(mu0)
        muPrev = mean(X,2);
else
        muPrev = mu0;
end
if nargin < 4 || isempty(sigma0)
        sigmaPrev = (X-muPrev)*(X-muPrev)'/n;
else
        sigmaPrev = sigma0;
end

%--------------------------------------------------------------
%%% CALCULATE ESTIMATE
%--------------------------------------------------------------

% Loop until converged
err = inf;
tmp = zeros(p,p,n);
while err > tol
    % Weights
    y = mdist(X,muPrev,sigmaPrev);
    % For some distributions, the MLE weights grow large near d=0, and this can result in instabilities
    % Use the following fudge factor that results in better numerical stability with negligible effect on the accuracy of the solution
    y(y<1e-10) = 1e-10;
    w = wFun(y);
    % Estimate mu
    muNext = sum(w.*X,2)/sum(w);
    % Estimate sigma
    for i = 1:n
        tmp(:,:,i) = (X(:,i)-muNext)*(X(:,i)-muNext)';
    end
    tmp2 = reshape(w,1,1,n) .* tmp;
    sigmaNext = mean(tmp2,3);
    % For now, error metric is Gaussian KL-divergence between previous and next estimates
    errPrev = err;
    if sum(isnan(sigmaNext(:)))
        error('error');
    end
    err = trace(sigmaNext\sigmaPrev) - p - log(det(sigmaPrev)/det(sigmaNext)) + (muNext-muPrev)'*(sigmaNext\(muNext-muPrev));

    % If error starts to grow, take a partial step
    xi0 = 0.9;
    xi = xi0;
    while err > errPrev    
        muNext = (1-xi)*muPrev + xi*muNext;
        sigmaNext = (1-xi)*sigmaPrev + xi*sigmaNext;
        err = trace(sigmaNext\sigmaPrev) - p - log(det(sigmaPrev)/det(sigmaNext)) + (muNext-muPrev)'*(sigmaNext\(muNext-muPrev));
        % Shrink xi in case another iteration is needed
        xi = xi0*xi;
        if xi < 0.001
            warning('Hit limit on xi');
        end
    end
    % Prepare for next iteration
    muPrev    = muNext;
    sigmaPrev = sigmaNext;
end
muHat    = muNext;
sigmaHat = sigmaNext;
omegaHat = sigmaHat/det(sigmaHat)^(1/p);
end
