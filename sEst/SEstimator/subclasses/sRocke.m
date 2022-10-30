% Copyright (C) 2020-2022 Justin A Fishbone
classdef sRocke < SEstimator % S-Estimator

    % Read-only properties
    properties (SetAccess = protected, GetAccess = public)
        isSq   = false;  % Is this an S-q estimator?
        maxTun = 1;      % Maximum tuning parameter value to use for enlargement logic
    end                  % (must result in: numerically computable result, proper S-estimator)
    
    % Public methods
    methods
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% Constructor
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function obj = sRocke(varargin)
        % obj = sRocke()
        % obj = sRocke(..., gamma, tol, b)  -- All optional inputs
            % Call superclass constructor
            obj = obj@SEstimator(varargin{:});
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% rhoTildFun(): Decreasing rho (rho-tilde) function
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function rho = rhoTildFun(obj, t, gam)
        % rho = rhoTildFun(obj, t)
        % rho = rhoTildFun(obj, t, gamma)
            if nargin < 3
                gam = obj.tun;
            end
            rho = (t-1)./(4*gam) .* (3 - ((t-1)./gam).^2) + 0.5;
            rho(t <= 1-gam) = 0;
            rho(t >= 1+gam) = 1;
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% wTildFun(): w-tilde function
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function w = wTildFun(obj, t)
            gam = obj.tun;
            w = 3/4/gam*(1-((t-1)./gam).^2); % Eq 6.42, Pg 212 from Maronna et al, 2019
            w(w<0) = 0;
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% extremesFun(): Estimator extremes
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [a, c, rhoA, rhoC] = extremesFun(obj, gam)
            if nargin < 2
                gam = obj.tun;
            end
            a    = 1 - gam;
            c    = 1 + gam;
            rhoA = 0;
            rhoC = 1;
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% distMedFun(): Return median of expected PDF
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function med = distMedFun(obj)
            if obj.isComplex
                k = 2;
            else
                k = 1;
            end
            med = chi2inv(0.5, k*p)/k; % Always assume Gaussian data for S-Rocke estimator
        end
    end
end