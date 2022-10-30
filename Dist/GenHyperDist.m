% Copyright (C) 2020-2022 Justin A Fishbone
classdef GenHyperDist
    % Static methods
    methods (Static)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% rand(): Generate n p-dimensional random vectors
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function X = rand(p, n, mu, sigma, lambda, psi, chi, isComplex)
        % X = GenHyperDist.rand(p, n, [],    [], lambda, psi, chi)
        % X = GenHyperDist.rand(p, n, mu, sigma, lambda, psi, chi)
        % X = GenHyperDist.rand(p, n, mu, sigma, lambda, psi, chi, isComplex)
        % 
        % INPUTS:
        % p         - Dimension
        % n         - Number of random vectors
        % mu        - (Optional) px1 location vector (default is zeros)
        % sigma     - (Optional) pxp scatter matrix (default is identity)
        % lambda    - Distribution parameter
        % psi       - Distribution parameter
        % chi       - Distribution parameter
        % isComplex - (Optional) Whether to generate complex-valued data (default is false)
        %
        % OUTPUTS:
        % X         - pxn n p-dimensional random vectors with elliptical generalized hyperbolic distribution
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
            if chi==0
                alpha = psi^(k*p/4+lambda/2)/pi^(k*p/2)/2^(k*p/2+lambda-1)/gamma(lambda);
            else
                alpha = psi^(k*p/4)/(2*pi)^(k*p/2)/chi^(lambda/2)/besselk(lambda,sqrt(chi*psi));
            end
            fr = @(r) sk * 2*pi^(k*p/2)/gamma(k*p/2)*alpha * (sk*r).^(k*p-1) .* besselk(lambda-k*p/2, sqrt(psi*(chi+(sk*r).^2)))./(chi+(sk*r).^2).^(k*p/4-lambda/2);
            % Other parameters
            TOL = 1e-10;
            R_REF_MAX_SURV = 2.^(1:0.5:13);
            N_REF_SAMPLES = 1e4;
            ABS_TOL = 1e-10;
            REL_TOL = 1e-10;
            % Call random elliptical function -- Note this cannot take nested function handles
            X = randEllip(fr, p, n, mu, sigma, R_REF_MAX_SURV, ABS_TOL, REL_TOL, TOL, N_REF_SAMPLES, isComplex);
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% mle(): Maximum likelihood estimator of location and scatter/shape
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [muHat, sigmaHat, omegaHat] = mle(X, mu0, sigma0, tol, lambda, psi, chi)
        % [muHat, sigmaHat, omegaHat] = GenHyperDist.mle(X, mu0, sigma0, tol, lambda, psi, chi)
        % Generalized hyperbolic MLE for location and scatter/shape of X
        %
        % INPUTS: 
        % X        - pxn n p-dimensional random vectors with assumed distribution
        % mu0      - (Optional) px1 initial guess of location
        % sigma0   - (Optional) pxp initial guess of scatter
        % tol      - (Optional) Convergence criteria for numerical solution of estimates
        % lambda   - Distribution parameter
        % psi      - Distribution parameter
        % chi      - Distribution parameter
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
            
            p = size(X,1);
            wFun = @(d) k * sqrt(psi./(chi+(k*d))).*besselk(lambda-p/2+1, sqrt(psi*(chi+(k*d))))./besselk(lambda-p/2, sqrt(psi*(chi+(k*d)))) - 2*(lambda-p/2)./(chi+(k*d));
            [muHat, sigmaHat, omegaHat] = solveMLE(X, wFun, mu0, sigma0, tol);
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% distMedFun(): Returns the median of the distribution of squared Mahalanobis distance
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function med = distMedFun(p, lambda, psi, chi, isComplex)
        % med = GenHyperDist.distMedFun(p, lambda, psi, chi, isComplex)
        % Returns the median of the distribution of squared Mahalanobis distance
        %
        % INPUTS: p, lambda, psi, chi - Distribution Parameters
        % isComplex - (Optional) Whether to generate complex-valued data (default is false)
        %
        
            if nargin < 5 || isempty(isComplex) || ~isComplex
                k = 1;
            else
                k = 2;
            end
            
            % Use persistent variables for efficiency since this may be 
            % called many times, and the value will rarely change
            persistent pp pMD2Median ll psps chichi
            if isempty(pp) || pp ~= p || ll ~= lambda || psps ~= psi || chichi ~= chi

                % PDF
                phi = @(d) besselk(lambda-k*p/2, sqrt(psi*(chi+d)))./(chi+d).^(k*p/4-lambda/2);
                if chi==0
                    alpha = psi^(k*p/4+lambda/2)/pi^(k*p/2)/2^(k*p/2+lambda-1)/gamma(lambda);
                else
                    alpha = psi^(k*p/4)/(2*pi)^(k*p/2)/chi^(lambda/2)/besselk(lambda,sqrt(chi*psi));
                end
                fd = @(d) k * pi^(k*p/2)/gamma(k*p/2)*alpha * (k*d).^(k*p/2-1).*phi(k*d);

                % Survey points to find approximate median
                R_REF_SURV = 2.^(-5:0.5:20);
                f_ref_surv = arrayfun( @(x) integral(fd, 0, x), R_REF_SURV);
                lowBInd = find(f_ref_surv < 0.5, 1, 'last');
                uppBInd = find(f_ref_surv > 0.5, 1, 'first');
                if isempty(lowBInd)
                    error('Unable to find lower bound');
                end
                if isempty(uppBInd)
                    error('Unable to find upper bound');
                end

                % Integrate PDF to get CDF, and find where CDF = 0.5
                pMD2Median = fzero( @(x) integral(fd, 0, x, 'AbsTol', 1e-12, 'RelTol', 1e-12) - 0.5, [R_REF_SURV(lowBInd), R_REF_SURV(uppBInd)]);
                pp = p;
                ll = lambda;
                psps = psi;
                chichi = chi;
            end
            med = pMD2Median;
        end
    end
end