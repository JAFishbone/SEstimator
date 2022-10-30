% Copyright (C) 2020-2022 Justin A Fishbone
classdef qPearsonVII < SEstimator % S-Estimator

    % Public properties
    properties
        s(1,1) double % Distribution parameter
        N(1,1) double % Distribution parameter
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
        function obj = qPearsonVII(s, N, varargin)
        % obj = qPearsonVII(s, N)
        % obj = qPearsonVII(..., tun0, tol, b, isComplex)  -- All optional inputs
            % Call superclass constructor
            obj = obj@SEstimator(varargin{:});
            
            % Set distribution parameters
            obj.s = s;
            obj.N = N;
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
            s = obj.s; N = obj.N; %#ok<*PROPLC>
            
            rho = -(gamma(N)/gamma(N-p/2)/gamma(p/2))^sq/s^(sq*p/2) * t.^(sp*sq) .* (1+t/s).^(-N*sq) .* (sp - N*t./(s+t));
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% wTildFun(): w-tilde function
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function w = wTildFun(obj, t)
            % For simplicity/clarity, pull properties into local variables
            p = obj.pEff; sp = p/2-1;
            sq = 1 - obj.tun;
            s = obj.s; N = obj.N;
            
            w = -(gamma(N)/gamma(N-p/2)/gamma(p/2))^sq/s^(sq*p/2) * t.^(sp*sq) .* (1+t/s).^(-N*sq) ...
                .* ( sq*(sp./sqrt(t) - N*sqrt(t)./(s+t)).^2 - N*s./(s+t).^2 );
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
            s = obj.s; N = obj.N;
            % Minimum location and corresponding value of rho-tilde
            a = s*( 2*N*sq*sp+N-2*sq*sp^2 - sqrt(4*N^2*sq*sp-4*N*sq*sp^2+N^2) ) ...
                / ( 2*sq*(sp^2-2*N*sp+N^2) );
            rhoA = obj.rhoTildFun(a,q);
            % Maximum location and corresponding value of rho-tilde
            c = s*( 2*N*sq*sp+N-2*sq*sp^2 + sqrt(4*N^2*sq*sp-4*N*sq*sp^2+N^2) ) ...
                / ( 2*sq*(sp^2-2*N*sp+N^2) );
            rhoC = obj.rhoTildFun(c,q);
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% distMedFun(): Return median of expected PDF
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function med = distMedFun(obj)
            med = PearsonVIIDist.distMedFun(obj.p, obj.s, obj.N, obj.isComplex);
        end
    end
end