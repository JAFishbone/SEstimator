% Copyright (C) 2020-2022 Justin A Fishbone
classdef PearsonIIDist
    % Static methods
    methods (Static)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% rand(): Generate n p-dimensional random vectors
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function X = rand(p, n, mu, sigma, m, isComplex)
        % X = PearsonIIDist.rand(p, n, [],    [], m)
        % X = PearsonIIDist.rand(p, n, mu, sigma, m)
        % X = PearsonIIDist.rand(p, n, mu, sigma, m, isComplex)
        % 
        % INPUTS:
        % p         - Dimension
        % n         - Number of random vectors
        % mu        - (Optional) px1 location vector (default is zeros)
        % sigma     - (Optional) pxp scatter matrix (default is identity)
        % m         - Distribution parameter
        % isComplex - (Optional) Whether to generate complex-valued data (default is false)
        %
        % OUTPUTS:
        % X         - pxn n p-dimensional random vectors with elliptical Pearson type II distribution
        %
            if nargin < 3 || isempty(mu)
                mu = zeros(p,1);
            end
            if nargin < 4 || isempty(sigma)
                sigma = eye(p);
            end
            if nargin < 6 || isempty(isComplex)
                isComplex = false;
            end
            
            % PDF Function
            if isComplex
                k = 2;
            else
                k = 1;
            end
            sk = sqrt(k);
            fr = @(r) sk * 2*gamma(k*p/2+m+1)/gamma(m+1)/gamma(k*p/2) * (sk*r).^(k*p-1) .* (1-(sk*r).^2).^m;
            % Other parameters
            R_REF_MAX_SURV = 1;
            % Call random elliptical function -- Note this cannot take nested function handles
            X = randEllip(fr, p, n, mu, sigma, R_REF_MAX_SURV, [], [], [], [], isComplex);
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% mle(): Maximum likelihood estimator of location and scatter/shape
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [muHat, sigmaHat, omegaHat] = mle(X, mu0, sigma0, tol, m)
        % [muHat, sigmaHat, omegaHat] = PearsonIIDist.mle(X, mu0, sigma0, tol, m)
        % Pearson type II MLE for location and scatter/shape of X
        %
        % INPUTS: 
        % X        - pxn n p-dimensional random vectors with assumed distribution
        % mu0      - (Optional) px1 initial guess of location
        % sigma0   - (Optional) pxp initial guess of scatter
        % tol      - (Optional) Convergence criteria for numerical solution of estimates
        % m        - Distribution parameter
        %
        % OUTPUTS:
        % muHat    - px1 location matrix estimate
        % sigmaHat - pxp scatter matrix estimate
        % omegaHat - pxp shape matrix estimate
        % 
        
            % Real or complex?
            if isreal(X)
                k = 1;
            else
                k = 2;
            end
            
            wFun = @(d) k * 2*m./(1-(k*d));
            [muHat, sigmaHat, omegaHat] = solveMLE(X, wFun, mu0, sigma0, tol);
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% distMedFun(): Returns the median of the distribution of squared Mahalanobis distance
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function med = distMedFun(p, m, isComplex)
        % med = PearsonIIDist.distMedFun(p, m, isComplex)
        % Returns the median of the distribution of squared Mahalanobis distance
        %
        % INPUTS: p, m - Distribution Parameters
        % isComplex - (Optional) Whether to generate complex-valued data (default is false)
        %
   
            if nargin < 3 || isempty(isComplex) || ~isComplex
                k = 1;
            else
                k = 2;
            end
            
            % Use persistent variables for efficiency since this may be 
            % called many times, and the value will rarely change
            persistent pp pMD2Median mm
            if isempty(pp) || pp ~= p || mm ~= m
                fd = @(d) k * gamma(k*p/2+m+1)/gamma(k*p/2)/gamma(m+1) * (k*d).^(k*p/2-1).*(1-(k*d)).^m;
                pMD2Median = fzero( @(x) integral(fd, 0, x) - 0.5, [1e-3, 1/k]);
                pp = p;
                mm = m;
            end
            med = pMD2Median;
        end
    end
end