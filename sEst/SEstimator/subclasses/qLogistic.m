% Copyright (C) 2020-2022 Justin A Fishbone
classdef qLogistic < SEstimator % S-Estimator

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
        function obj = qLogistic(varargin)
        % obj = qLogistic()
        % obj = qLogistic(..., tun0, tol, b, isComplex)  -- All optional inputs
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
            if nargin < 3
                q = obj.tun; %#ok<*PROPLC>
            end
            sq = 1 - q;
            sp = obj.pEff/2-1;
            rho = -t.^(sq*sp).*exp(-t*sq)./(1+exp(-t)).^(2*sq) .* (sp + t.*(exp(-t)-1)./(exp(-t)+1));
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% wTildFun(): w-tilde function
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function w = wTildFun(obj, t, q)
        % w = wTildFun(obj, t)
        % w = wTildFun(obj, t, q)
            % For simplicity/clarity, pull properties into local variables
            % For simplicity/clarity, pull properties into local variables
            if nargin < 3
                q = obj.tun;
            end
            sq = 1 - q;
            sp = obj.pEff/2-1;
            % Intermediate value
            pp = exp(-t)./(1+exp(-t));

            % Weight function. Drop scaling constants for simplicity (and speed)
            w = -t.^(sq*sp).*exp(-t*sq)./(1+exp(-t)).^(2*sq) ...
                .* ( sq*sp^2./t + (2*sq*sp+1)*(2*pp-1) -q*t.*(2*pp-1).^2 ...
                      + t.*(8*pp.^2 - 8*pp.^2 - 6*pp./(1+exp(-t)) + 1) );
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% extremesFun(): Estimator extremes
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [a, c, rhoA, rhoC] = extremesFun(obj, q)
            % For simplicity/clarity, pull properties into local variables
            p = obj.pEff;
            if nargin < 2
                q = obj.tun;
            end
            
            % Use persistent variables for efficiency since this may be 
            % called many times, and the value will rarely change
            persistent pp qq aa cc rhoAA rhoCC;
            if isempty(aa) || isnan(aa) || isempty(cc) || isnan(cc) || pp ~= p || qq ~= q

                % Do a quick survey of equally-spaced points to narrow-in on zeros
                ySurvey = 0.001:100;
                isPositive = obj.wTildFun(ySurvey,q) > 0;
                % Indices into ySurvey: 
                %     ind1, ind2 should straddle the first zero, 
                %     and ind3, ind4 should straddle the second (last) zero
                ind2 = find(isPositive, 1, 'first');
                ind1 = find(~isPositive(1:ind2), 1, 'last');
                ind3 = (ind2-1) + find(isPositive(ind2:end), 1, 'last');
                ind4 = (ind2-1) + find(~isPositive(ind2:end), 1, 'first');
                if ~isempty(ind4)
                    yEnd = ySurvey(ind4);
                else
                    yTest = ySurvey(end);
                    yVal  = 1;
                    while yVal > 0
                        yTest = 5*yTest;
                        yVal = obj.wTildFun(yTest,q);
                    end
                    yEnd = yTest;
                end

                % Search for the two zeros
                aa = fzero( @(y) obj.wTildFun(y,q), [ySurvey(ind1), ySurvey(ind2)]);
                cc = fzero( @(y) obj.wTildFun(y,q), [ySurvey(ind3),          yEnd] );
                if aa==cc
                    error('Extrema expected to differ');
                end
                pp = p;
                qq = q;
                rhoAA = obj.rhoTildFun(aa,q);
                rhoCC = obj.rhoTildFun(cc,q);
            end
            a    = aa;
            c    = cc;
            rhoA = rhoAA;
            rhoC = rhoCC;
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% distMedFun(): Return median of expected PDF
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function med = distMedFun(obj)
            med = LogisticDist.distMedFun(obj.p, obj.isComplex);
        end
    end
end