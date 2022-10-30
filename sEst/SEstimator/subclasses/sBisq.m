% Copyright (C) 2020-2022 Justin A Fishbone
classdef sBisq < SEstimator % S-Estimator

    % Read-only properties
    properties (SetAccess = protected, GetAccess = public)
        isSq   = false;  % Is this an S-q estimator?
        maxTun = 1;      % Maximum tuning parameter value to use for enlargement logic
    end                  % (not applicable for bisquare, but set to 1 anyway)
    
    % Public methods
    methods
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% Constructor
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function obj = sBisq(varargin)
        % obj = sBisq()
        % obj = sBisq(..., tun0, tol, b)  -- All optional inputs, tuning is N/A to set to 1
            
            % Tuning parameter is N/A so set it to 1
            if nargin > 0
                varargin{1} = 1;
            end
            % Call superclass constructor
            obj = obj@SEstimator(varargin{:});
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% rhoTildFun(): Decreasing rho (rho-tilde) function
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function rho = rhoTildFun(~, t, ~)
            rho = min(1, 1-(1-t).^3);
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% wTildFun(): w-tilde function
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function w = wTildFun(~, t)
            w = 3*(1-t).^2.*(t<=1);
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% extremesFun(): Estimator extremes
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [a, c, rhoA, rhoC] = extremesFun(~, ~)
            a    = 0;
            c    = 1;
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
            med = chi2inv(0.5, k*p)/k; % Always assume Gaussian data for S-Bisquare estimator
        end
    end
end