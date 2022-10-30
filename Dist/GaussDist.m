% Copyright (C) 2020-2022 Justin A Fishbone
classdef GaussDist
    % Static methods
    methods (Static)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% rand(): Generate n p-dimensional random vectors
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function X = rand(p, n, mu, sigma, isComplex)
        % X = GaussDist.rand(p, n)
        % X = GaussDist.rand(p, n, mu, sigma)
        % X = GaussDist.rand(p, n, mu, sigma, isComplex)
        % 
        % INPUTS:
        % p         - Dimension
        % n         - Number of random vectors
        % mu        - (Optional) px1 location vector (default is zeros)
        % sigma     - (Optional) pxp scatter matrix (default is identity)
        % isComplex - (Optional) Whether to generate complex-valued data (default is false)
        %
        % OUTPUTS:
        % X         - pxn n p-dimensional random vectors with assumed distribution
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
            
            if isComplex
                X = mu + cholcov(sigma)'*complex(randn(p,n),randn(p,n))/sqrt(2);
            else
                X = mu + cholcov(sigma)'*randn(p,n);
            end
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% mle(): Maximum likelihood estimator of location and scatter/shape
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [muHat, sigmaHat, omegaHat] = mle(X)
        % [muHat, sigmaHat, omegaHat] = GaussDist.mle(X)
        % Gaussian MLE for location and scatter/shape of X
        %
        % INPUTS: 
        % X        - pxn n p-dimensional random vectors with assumed distribution
        %
        % OUTPUTS:
        % muHat    - px1 location matrix estimate
        % sigmaHat - pxp scatter matrix estimate
        % omegaHat - pxp shape matrix estimate
        % 
            p = size(X,1);
            n = size(X,2);
            muHat    = mean(X,2);
            sigmaHat = (X-muHat)*(X-muHat)'/n;
            omegaHat = sigmaHat/det(sigmaHat)^(1/p);
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% distMedFun(): Returns the median of the distribution of squared Mahalanobis distance
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function med = distMedFun(p, isComplex)
        % med = GaussDist.distMedFun(p, isComplex)
        % Returns the median of the distribution of squared Mahalanobis distance
        %
        % INPUTS: 
        % p         - Dimension
        % isComplex - (Optional) Whether to generate complex-valued data (default is false)
        %
        
            if nargin < 2 || isempty(isComplex) || ~isComplex
                k = 1;
            else
                k = 2;
            end
            
            med = chi2inv(0.5, k*p)/k;
        end
    end
end