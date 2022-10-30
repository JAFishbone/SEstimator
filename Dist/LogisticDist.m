% Copyright (C) 2020-2022 Justin A Fishbone
classdef LogisticDist
    % Static methods
    methods (Static)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% rand(): Generate n p-dimensional random vectors
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function X = rand(p, n, mu, sigma, isComplex)
        % X = LogisticDist.rand(p, n)
        % X = LogisticDist.rand(p, n, mu, sigma)
        % X = LogisticDist.rand(p, n, mu, sigma, isComplex)
        % 
        % INPUTS:
        % p         - Dimension
        % n         - Number of random vectors
        % mu        - (Optional) px1 location vector (default is zeros)
        % sigma     - (Optional) pxp scatter matrix (default is identity)
        % isComplex - (Optional) Whether to generate complex-valued data (default is false)
        %
        % OUTPUTS:
        % X         - pxn n p-dimensional random vectors with elliptical logistic distribution
        %
            if nargin < 3 || isempty(mu)
                mu = zeros(p,1);
            end
            if nargin < 4 || isempty(sigma)
                sigma = eye(p);
            end
            if nargin < 5 || isempty(isComplex)
                isComplex = false;
            end
            
            % PDF Function
            if isComplex
                k = 2;
            else
                k = 1;
            end
            sk = sqrt(k);
            phi = @(y) exp(-y)./(1+exp(-y)).^2;
            alpha = gamma(k*p/2)/pi^(k*p/2)/integral( @(y) y.^(k*p/2-1).*phi(y), 0, inf);
            fr = @(r) sk * 2*pi^(k*p/2)/gamma(k*p/2)*alpha * (sk*r).^(k*p-1).*exp(-(sk*r).^2)./(1+exp(-(sk*r).^2)).^2;
            % Other parameters
            TOL = 1e-12;
            R_REF_MAX_SURV = 11;
            N_REF_SAMPLES = 1e4;
            ABS_TOL = 1e-10;
            REL_TOL = 1e-6;
            % Call random elliptical function -- Note this cannot take nested function handles
            X = randEllip(fr, p, n, mu, sigma, R_REF_MAX_SURV, ABS_TOL, REL_TOL, TOL, N_REF_SAMPLES, isComplex);
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% mle(): Maximum likelihood estimator of location and scatter/shape
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [muHat, sigmaHat, omegaHat] = mle(X, mu0, sigma0, tol)
        % [muHat, sigmaHat, omegaHat] = LogisticDist.mle(X, mu0, sigma0, tol)
        % Logistic MLE for location and scatter/shape of X
        %
        % INPUTS: 
        % X        - pxn n p-dimensional random vectors with assumed distribution
        % mu0      - (Optional) px1 initial guess of location
        % sigma0   - (Optional) pxp initial guess of scatter
        % tol      - (Optional) Convergence criteria for numerical solution of estimates
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
            
            wFun = @(d) k * 2-4*exp(-(k*d))./(1+exp(-(k*d)));
            [muHat, sigmaHat, omegaHat] = solveMLE(X, wFun, mu0, sigma0, tol);
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% distMedFun(): Returns the median of the distribution of squared Mahalanobis distance
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function med = distMedFun(p, isComplex)
        % med = LogisticDist.distMedFun(p, isComplex)
        % Returns the median of the distribution of squared Mahalanobis distance
        %
        % INPUTS: p - Dimension
        % isComplex - (Optional) Whether to generate complex-valued data (default is false)
        %
        
            if nargin < 2 || isempty(isComplex) || ~isComplex
                k = 1;
            else
                k = 2;
            end
            % Use persistent variables for efficiency since this may be 
            % called many times, and the value will rarely change
            persistent pp pMD2Median;
            if isempty(pp) || pp ~= p
                a = 1/integral(@(y) y.^(k*p/2-1).*exp(-y)./(1+exp(-y)).^2, 0, inf);
                fd = @(d) k * a*(k*d).^(k*p/2-1).*exp(-(k*d))./(1+exp(-(k*d))).^2;
                pMD2Median = fzero( @(x) integral(fd, 0, x) - 0.5, [1, 100]);
                pp = p;
            end
            med = pMD2Median;
        end
    end
end