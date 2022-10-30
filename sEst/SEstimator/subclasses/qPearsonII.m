% Copyright (C) 2020-2022 Justin A Fishbone
classdef qPearsonII < SEstimator % S-Estimator

    % Public properties
    properties
        m(1,1) double % Distribution parameter
    end
    
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
        function obj = qPearsonII(m, varargin)
        % obj = qPearsonII(m)
        % obj = qPearsonII(..., tun0, tol, b, isComplex)  -- All optional inputs
            % Call superclass constructor
            obj = obj@SEstimator(varargin{:});
            
            % Set distribution parameters
            obj.m = m;
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% rhoTildFun(): Decreasing rho (rho-tilde) function
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function rho = rhoTildFun(obj, t, q)
        % rho = rhoTildFun(obj, t)
        % rho = rhoTildFun(obj, t, q)
            % For simplicity/clarity, pull properties into local variables
            p = obj.pEff; sp = p/2-1;
            if nargin < 3
                q = obj.tun;
            end
            sq = 1 - q;
            m = obj.m; %#ok<*PROPLC>
            
            rho = -(gamma(p/2+m+1)/gamma(p/2)/gamma(m+1))^sq * t.^(sq*sp).*(1-t).^(m*sq) .* (sp - m*t./(1-t));
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% wTildFun(): w-tilde function
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function w = wTildFun(obj, t)
            % For simplicity/clarity, pull properties into local variables
            p = obj.pEff; sp = p/2-1;
            sq = 1 - obj.tun;
            m = obj.m;
            
            w = -(gamma(p/2+m+1)/gamma(p/2)/gamma(m+1))^sq * t.^(sq*sp).*(1-t).^(m*sq) ...
                .* ( sq*sp^2./t - m*(2*sq*sp+1)./(1-t) + t./(1-t).^2.*(m^2*sq-m) );
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% extremesFun(): Estimator extremes
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [a, c, rhoA, rhoC] = extremesFun(obj, q)
            % For simplicity/clarity, pull properties into local variables
            if nargin < 2
                q = obj.tun;
            end
            sp = obj.pEff/2-1; sq = 1 - q;
            m = obj.m;
            
            % Check Conditions
            if q >= 1 - 1/m && q ~= 1
                error('Invalid q');
            end
            
            % Minimum location and corresponding value of rho-tilde
            a = ( 2*sq*sp^2+m*(2*sq*sp+1) - sqrt(m^2*(4*sq*sp+1)+4*m*sq*sp^2) ) ...
                / (2*(sq*sp^2+m*(2*sq*sp+m*sq)));
            rhoA = obj.rhoTildFun(a,q);
            % Maximum location and corresponding value of rho-tilde
            c = ( 2*sq*sp^2+m*(2*sq*sp+1) + sqrt(m^2*(4*sq*sp+1)+4*m*sq*sp^2) ) ...
                / (2*(sq*sp^2+m*(2*sq*sp+m*sq)));
            rhoC = obj.rhoTildFun(c,q);
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% distMedFun(): Return median of expected PDF
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function med = distMedFun(obj)
            med = PearsonIIDist.distMedFun(obj.p, obj.m, obj.isComplex);
        end
    end
end