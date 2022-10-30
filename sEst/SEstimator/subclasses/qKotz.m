% Copyright (C) 2020-2022 Justin A Fishbone
classdef qKotz < SEstimator % S-Estimator

    % Public properties
    properties
        s(1,1) double % Distribution parameter
        N(1,1) double % Distribution parameter
        r(1,1) double % Distribution parameter
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
        function obj = qKotz(s, N, r, varargin)
        % obj = qKotz(s,N,r)
        % obj = qKotz(..., tun0, tol, b, isComplex)  -- All optional inputs
            % Call superclass constructor
            obj = obj@SEstimator(varargin{:});
            
            % Set distribution parameters
            obj.s = s;
            obj.N = N;
            obj.r = r;
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% rhoTildFun(): Decreasing rho (rho-tilde) function
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function rho = rhoTildFun(obj, t, q)
        % rho = rhoTildFun(obj, t)
        % rho = rhoTildFun(obj, t, q)
            % For simplicity/clarity, pull properties into local variables
            p = obj.pEff;
            if nargin < 3
                q = obj.tun;
            end
            sq = 1 - q;
            s = obj.s; N = obj.N; r = obj.r; %#ok<*PROPLC>
            rho = -(s/gamma((2*N+p)/2/s)*r^((2*N+p)/2/s))^sq ...
                  * t.^(sq*(p/2-1+N)).*exp(-sq*r*t.^s) ...
                  .* ((p/2-1)+N-r*s*t.^s);
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% wTildFun(): w-tilde function
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function w = wTildFun(obj, t)
            % For simplicity/clarity, pull properties into local variables
            sp = obj.pEff/2 - 1; sq = 1 - obj.tun;
            s = obj.s; N = obj.N; r = obj.r; %#ok<*PROP>
            % Weight function
            w = -(s/gamma((1+N+sp)/s)*r^((1+N+sp)/s))^sq * t.^(sq*(sp+N)-1).*exp(-sq*r*t.^s) ...
                .* (sq*sp^2 + 2*sp*sq*N + sq*N^2 ...
                    + t.^s .* (sq*r^2*s^2*t.^s - r*s^2 - 2*r*s*sq*N - 2*r*s*sq*sp) );
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
            s = obj.s; N = obj.N; r = obj.r;

            % Check Conditions
            if N < -sp && q <= 1 + s/(4*(sp+N))
                error('Invalid q');
            end

            % Evaluate extremes
            ysa = (r*s^2 + 2*r*s*sq*N + 2*r*s*sq*sp - r*s*sqrt(s^2 + 4*s*sq*N + 4*s*sp*sq)) / (2*sq*r^2*s^2);
            ysc = (r*s^2 + 2*r*s*sq*N + 2*r*s*sq*sp + r*s*sqrt(s^2 + 4*s*sq*N + 4*s*sp*sq)) / (2*sq*r^2*s^2);
            a = ysa^(1/s);
            c = ysc^(1/s);
            rhoA = obj.rhoTildFun(a,q);
            rhoC = obj.rhoTildFun(c,q);
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% distMedFun(): Return median of expected PDF
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function med = distMedFun(obj)
            med = KotzDist.distMedFun(obj.p, obj.s, obj.N, obj.r, obj.isComplex);
        end
    end
end