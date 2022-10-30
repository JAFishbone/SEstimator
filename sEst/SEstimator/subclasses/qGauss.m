% Copyright (C) 2020-2022 Justin A Fishbone
classdef qGauss < SEstimator % S-Estimator

    % Read-only properties
    properties (SetAccess = protected, GetAccess = public)
        isSq   = true;   % Is this an S-q estimator?
        maxTun = 0.9999; % Maximum tuning parameter value to use for enlargement logic
    end                  % (must result in: numerically computable result, proper S-estimator)
    
    % Public methods
    methods
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% Constructor
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function obj = qGauss(varargin)
        % obj = qGauss()
        % obj = qGauss(q, tol, b, isComplex)  -- All optional inputs
            % Call superclass constructor
            obj = obj@SEstimator(varargin{:});
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% rhoTildFun(): Decreasing rho (rho-tilde) function
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function rho = rhoTildFun(obj, t, q)
        % rho = rhoTildFun(obj, t)
        % rho = rhoTildFun(obj, t, q)
            % For simplicity/clarity, pull properties into local variables
            p = obj.pEff; sp = obj.pEff/2 - 1; 
            if nargin > 2
                sq = 1 - q;
            else
                sq = 1 - obj.tun;
            end
            rho = (t/2-sp) .* t.^(sp*sq).*exp(-t/2*sq)/2^(p/2*sq)/gamma(p/2)^sq;
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% wTildFun(): w-tilde function
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function w = wTildFun(obj, t)
            % For simplicity/clarity, pull properties into local variables
            p = obj.pEff; sp = obj.pEff/2 - 1; sq = 1 - obj.tun;
            w = (1/2+sq*(sp-t/4-sp^2./t)) .* t.^(sp*sq).*exp(-t/2*sq)/2^(p/2*sq)/gamma(p/2)^sq;
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% extremesFun(): Estimator extremes
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [a, c, rhoA, rhoC] = extremesFun(obj, q)
            % For simplicity/clarity, pull properties into local variables
            sp = obj.pEff/2 - 1;
            if nargin < 2
                q = obj.tun;
            end
            sq = 1 - q;
            % Minimum location and corresponding value of rho-tilde
            a = (1+2*sp*sq-sqrt(1+4*sp*sq))/sq;
            rhoA = obj.rhoTildFun(a,q);
            % Maximum location and corresponding value of rho-tilde
            c = (1+2*sp*sq+sqrt(1+4*sp*sq))/sq;
            rhoC = obj.rhoTildFun(c,q);
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% distMedFun(): Return median of expected PDF
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function med = distMedFun(obj)
            med = GaussDist.distMedFun(obj.p, obj.isComplex);
        end
    end
end