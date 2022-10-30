% Copyright (C) 2020-2022 Justin A Fishbone
classdef PearsonVIIDist
    % Static methods
    methods (Static)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% rand(): Generate n p-dimensional random vectors
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function X = rand(p, n, mu, sigma, s, N, isComplex)
        % X = PearsonVIIDist.rand(p, n, [],    [], s, N)
        % X = PearsonVIIDist.rand(p, n, mu, sigma, s, N)
        % X = PearsonVIIDist.rand(p, n, mu, sigma, s, N, isComplex)
        % 
        % INPUTS:
        % p         - Dimension
        % n         - Number of random vectors
        % mu        - (Optional) px1 location vector (default is zeros)
        % sigma     - (Optional) pxp scatter matrix (default is identity)
        % s         - Distribution parameter
        % N         - Distribution parameter
        % isComplex - (Optional) Whether to generate complex-valued data (default is false)
        %
        % OUTPUTS:
        % X         - pxn n p-dimensional random vectors with elliptical Pearson type VII distribution
        %
            if nargin < 3 || isempty(mu)
                mu = zeros(p,1);
            end
            if nargin < 4 || isempty(sigma)
                sigma = eye(p);
            end
            if nargin < 7 || isempty(isComplex)
                isComplex = false;
            end
            
            % PDF Function
            if isComplex
                k = 2;
            else
                k = 1;
            end
            sk = sqrt(k);
            fr = @(r) sk * 2*gamma(N)/gamma(N-k*p/2)/s^(k*p/2)/gamma(k*p/2) * (sk*r).^(k*p-1) .* (1 + (sk*r).^2/s).^(-N);
            % Other parameters
            TOL = 1e-9;
            R_REF_MAX_SURV = 2.^(0:0.5:40);
            % Call random elliptical function -- Note this cannot take nested function handles
            X = randEllip(fr, p, n, mu, sigma, R_REF_MAX_SURV, [], [], TOL, [], isComplex);
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% mle(): Maximum likelihood estimator of location and scatter/shape
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [muHat, sigmaHat, omegaHat] = mle(X, mu0, sigma0, tol, s, N)
        % [muHat, sigmaHat, omegaHat] = PearsonVIIDist.mle(X, mu0, sigma0, tol, s, N)
        % Pearson type VII MLE for location and scatter/shape of X
        %
        % INPUTS: 
        % X        - pxn n p-dimensional random vectors with assumed distribution
        % mu0      - (Optional) px1 initial guess of location
        % sigma0   - (Optional) pxp initial guess of scatter
        % tol      - (Optional) Convergence criteria for numerical solution of estimates
        % s        - Distribution parameter
        % N        - Distribution parameter
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
            
            wFun = @(d) k * 2*N./(s+(k*d));
            [muHat, sigmaHat, omegaHat] = solveMLE(X, wFun, mu0, sigma0, tol);
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% distMedFun(): Returns the median of the distribution of squared Mahalanobis distance
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function med = distMedFun(p, s, N, isComplex)
        % med = PearsonVIIDist.distMedFun(p, s, N, isComplex)
        % Returns the median of the distribution of squared Mahalanobis distance
        %
        % INPUTS: p, s, N - Distribution Parameters
        % isComplex - (Optional) Whether to generate complex-valued data (default is false)
        %
        
            if nargin < 4 || isempty(isComplex) || ~isComplex
                k = 1;
            else
                k = 2;
            end

            % Use persistent variables for efficiency since this may be 
            % called many times, and the value will rarely change
            persistent pp pMD2Median ss NN
            if isempty(pp) || pp ~= p || ss ~= s || NN ~= N
                fd = @(d) k * gamma(N)/gamma(N-k*p/2)/gamma(k*p/2)/s^(k*p/2) * (k*d).^(k*p/2-1) .* (1+(k*d)/s).^(-N);
                pMD2Median = fzero( @(x) integral(fd, 0, x) - 0.5, [1e-3, 1e3]);
                pp = p;
                ss = s;
                NN = N;
            end
            med = pMD2Median;
        end
    end
end