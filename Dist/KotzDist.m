% Copyright (C) 2020-2022 Justin A Fishbone
classdef KotzDist
    % Static methods
    methods (Static)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% rand(): Generate n p-dimensional random vectors
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function X = rand(p, n, mu, sigma, s, N, r, isComplex)
        % X = KotzDist.rand(p, n, [],    [], s, N, r)
        % X = KotzDist.rand(p, n, mu, sigma, s, N, r)
        % X = KotzDist.rand(p, n, mu, sigma, s, N, r, isComplex)
        % 
        % INPUTS:
        % p         - Dimension
        % n         - Number of random vectors
        % mu        - (Optional) px1 location vector (default is zeros)
        % sigma     - (Optional) pxp scatter matrix (default is identity)
        % s         - Distribution parameter
        % N         - Distribution parameter
        % r         - Distribution parameter
        % isComplex - (Optional) Whether to generate complex-valued data (default is false)
        %
        % OUTPUTS:
        % X         - pxn n p-dimensional random vectors with elliptical Kotz type distribution
        %
            if nargin < 3 || isempty(mu)
                mu = zeros(p,1);
            end
            if nargin < 4 || isempty(sigma)
                sigma = eye(p);
            end
            if nargin < 8 || isempty(isComplex)
                isComplex = false;
            end

            % PDF Function
            if isComplex
                k = 2;
            else
                k = 1;
            end
            sk = sqrt(k);
            alpha = s*gamma(k*p/2)/pi^(k*p/2)/gamma((2*N+k*p)/2/s)*r^((2*N+k*p)/2/s);
            fr = @(y) sk * 2*pi^(k*p/2)/gamma(k*p/2)*alpha * (sk*y).^(k*p-1) .* (sk*y).^(2*N) .* exp(-r*(sk*y).^(2*s));

            % Other parameters
            TOL = 1e-10;
            R_REF_MAX_SURV = 2.^(0:45);
            N_REF_SAMPLES = 1e4;
            ABS_TOL = 1e-12;
            REL_TOL = 1e-12;
            % Call random elliptical function -- Note this cannot take nested function handles
            X = randEllip(fr, p, n, mu, sigma, R_REF_MAX_SURV, ABS_TOL, REL_TOL, TOL, N_REF_SAMPLES, isComplex);
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% mle(): Maximum likelihood estimator of location and scatter/shape
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [muHat, sigmaHat, omegaHat] = mle(X, mu0, sigma0, tol, s, N, r)
        % [muHat, sigmaHat, omegaHat] = KotzDist.mle(X, mu0, sigma0, tol, s, N, r)
        % Kotz type MLE for location and scatter/shape of X
        %
        % INPUTS: 
        % X        - pxn n p-dimensional random vectors with assumed distribution
        % mu0      - (Optional) px1 initial guess of location
        % sigma0   - (Optional) pxp initial guess of scatter
        % tol      - (Optional) Convergence criteria for numerical solution of estimates
        % s        - Distribution parameter
        % N        - Distribution parameter
        % r        - Distribution parameter
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
            
            wFun = @(d) k * 2*r*s*(k*d).^(s-1) - 2*N./(k*d);
            [muHat, sigmaHat, omegaHat] = solveMLE(X, wFun, mu0, sigma0, tol);
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% distMedFun(): Returns the median of the distribution of squared Mahalanobis distance
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function med = distMedFun(p, s, N, r, isComplex)
        % med = KotzDist.distMedFun(p, s, N, r, isComplex)
        % Returns the median of the distribution of squared Mahalanobis distance
        %
        % INPUTS: p, s, N, r - Distribution Parameters
        % isComplex - (Optional) Whether to generate complex-valued data (default is false)
        %
        
            if nargin < 5 || isempty(isComplex) || ~isComplex
                k = 1;
            else
                k = 2;
            end
        
            % Use persistent variables for efficiency since this may be 
            % called many times, and the value will rarely change
            persistent pp pMD2Median ss NN rr
            if isempty(pp) || pp ~= p || ss ~= s || NN ~= N || rr ~= r
                warning('off', 'MATLAB:integral:NonFiniteValue');
                % PDF
                phi = @(d) d.^N .* exp(-r*d.^s);
                alpha = s*gamma(k*p/2)/pi^(k*p/2)/gamma((2*N+k*p)/2/s)*r^((2*N+k*p)/2/s);
                fd = @(d) k * pi^(k*p/2)/gamma(k*p/2)*alpha * (k*d).^(k*p/2-1).*phi(k*d);

                % Survey points to find approximate median
                R_REF_SURV = 2.^(-5:0.5:80);
                f_ref_surv = arrayfun( @(x) integral(fd, 0, x), R_REF_SURV);
                uppBInd = find(f_ref_surv > 0.5, 1, 'first');
                if isempty(uppBInd)
                    error('Unable to find upper bound');
                end

                % Integrate PDF to get CDF, and find where CDF = 0.5
                pMD2Median = fzero( @(x) integral(fd, 0, x, 'AbsTol', 1e-12, 'RelTol', 1e-12) - 0.5, [R_REF_SURV(uppBInd-1), R_REF_SURV(uppBInd)]);
                pp = p;
                ss = s;
                NN = N;
                rr = r;
                warning('on', 'MATLAB:integral:NonFiniteValue');
            end
            med = pMD2Median;
        end
    end
end