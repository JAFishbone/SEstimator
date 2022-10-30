% Copyright (C) 2020-2022 Justin A Fishbone
classdef qGenHyper < SEstimator % S-Estimator

    % Public properties
    properties
        lambda(1,1) double % Distribution parameter
        psi(1,1)    double % Distribution parameter
        chi(1,1)    double % Distribution parameter
    end
    
    % Read-only properties
    properties (SetAccess = protected, GetAccess = public)
        isSq   = true;   % Is this an S-q estimator?
        maxTun = 0.998;  % Above ~0.998 rho and w functions lack numerical sufficient
    end                  % precision due to modified Bessel function of 2nd kind
    
    % Public methods
    methods
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% Constructor
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function obj = qGenHyper(lambda, psi, chi, varargin)
        % obj = qGenHyper(lambda, psi, chi)
        % obj = XXX(..., tun0, tol, b, isComplex)  -- All optional inputs
            % Call superclass constructor
            obj = obj@SEstimator(varargin{:});
            
            % Set distribution parameters
            obj.lambda = lambda;
            obj.psi = psi;
            obj.chi = chi;
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% rhoTildFun(): Decreasing rho (rho-tilde) function
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function rho = rhoTildFun(obj, t, q)
        % rho = rhoTildFun(obj, t)
        % rho = rhoTildFun(obj, t, q)
            % For simplicity/clarity, pull properties into local variables
            if nargin < 3
                q = obj.tun;
            end
            sq = 1 - q;
            p = obj.pEff;
            lambda = obj.lambda; psi = obj.psi; chi = obj.chi; %#ok<*PROPLC>
            % Scale factor
            if chi == 0
                beta = psi^(sq*(p/4+lambda/2))/2^(sq*(p/2+lambda-1))/(gamma(lambda)*gamma(p/2))^sq;
            else
                beta = psi^(sq*p/4)/2^(sq*p/2)/chi^(sq*lambda/2)/(gamma(p/2)*besselk(lambda,sqrt(chi*psi)))^sq;
            end
            % Decreasing rho function
            rho = -beta * t.^(sq*(p/2-1)) .* besselk(lambda-p/2, sqrt(psi*(chi+t))).^sq ./ (chi+t).^(sq*(p/4-lambda/2)) ...
                      .* (p/2-1 + t.* ( (lambda-p/2)./(chi+t) ...
                                        - sqrt(psi)*besselk(lambda-p/2+1, sqrt(psi*(chi+t)))/2./sqrt(chi+t)./besselk(lambda-p/2, sqrt(psi*(chi+t)))));
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% wTildFun(): w-tilde function
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function w = wTildFun(obj, t, q)
        % w = wTildFun(obj, t)
        % w = wTildFun(obj, t, q)
            % For simplicity/clarity, pull properties into local variables
            if nargin < 3
                q = obj.tun;
            end
            sq = 1 - q;
            % For simplicity/clarity, pull properties into local variables
            p = obj.pEff;
            lambda = obj.lambda; psi = obj.psi; chi = obj.chi;
            
            % Simplified constants
            sp = p/2 - 1;
            sl = lambda - p/2;
            % Intermediate values
            cpt = chi + t;
            bl  = besselk(sl,   sqrt(psi*(cpt)));
            bl1 = besselk(sl+1, sqrt(psi*(cpt)));
            % Scale factor
            if chi == 0
                beta = psi^(sq*(p/4+lambda/2))/2^(sq*(p/2+lambda-1))/(gamma(lambda)*gamma(p/2))^sq;
            else
                beta = psi^(sq*p/4)/2^(sq*p/2)/chi^(sq*lambda/2)/(gamma(p/2)*besselk(lambda,sqrt(chi*psi)))^sq;
            end
            % Weight function
            w = -beta * t.^(sq*sp) .* bl.^sq .* cpt.^(sq*sl/2) ...
                .* ( sq*sp^2./t ...
                     + 1./cpt .* (sl*(2*sq*sp+1) + t .* ( sl*(sl-1)./cpt - q*sl^2./cpt + psi/4 ) ) ...
                     + sqrt(psi)*bl1./sqrt(cpt)./bl .* ( t.*(sl*(2*q-1)+1)/2./cpt - q*sqrt(psi).*t/4./sqrt(cpt).*bl1./bl - sq*sp-1/2 ) );
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% extremesFun(): Estimator extremes
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [a, c, rhoA, rhoC] = extremesFun(obj, q)
            % For simplicity/clarity, pull properties into local variables
            p = obj.pEff;
            lambda = obj.lambda; psi = obj.psi; chi = obj.chi;
            if nargin < 2
                q = obj.tun;
            end
            
            % Use persistent variables for efficiency since this may be 
            % called many times, and the value will rarely change
            persistent pp qq aa cc rhoAA rhoCC ll psps chichi
            if isempty(aa) || isnan(aa) || isempty(cc) || isnan(cc) || pp ~= p || qq ~= q || ll ~= lambda || psps ~= chi || chichi ~= chi
                % Do a quick survey of equally-spaced points to narrow-in on zeros
                ySurvey = 10.^(-10:0.025:10);
                samples = obj.wTildFun(ySurvey, q);
                failedSamples = isnan(samples) | isinf(samples);
                isPositive = samples > 0;
                % Indices into ySurvey: 
                %     ind1, ind2 should straddle the first zero, 
                %     and ind3, ind4 should straddle the second (last) zero
                ind2 = find(isPositive & ~failedSamples, 1, 'first');
                if isempty(ind2)
                    error('Unable to find positive weight sample value');
                end
                ind1 = find(~isPositive(1:ind2) & ~failedSamples(1:ind2), 1, 'last');
                ind4 = (ind2-1) + find(~isPositive(ind2:end) & ~failedSamples(ind2:end), 1, 'first');
                ind3 = (ind2-1) + find(isPositive(ind2:ind4) & ~failedSamples(ind2:ind4), 1, 'last');    
                if isempty(ind4)
                    error('Weight function not changing signs at top');
                end

                % Search for the two zeros
                if isempty(ind1)
                    % warning('Weight function not changing signs at bottom, a=0');
                    aa = 0;
                else
                    aa = fzero( @(y) obj.wTildFun(y,q), [ySurvey(ind1), ySurvey(ind2)]);
                end
                cc = fzero( @(y) obj.wTildFun(y,q), [ySurvey(ind3),ySurvey(ind4)] );
                if aa==cc
                    error('Extrema expected to differ');
                end
                pp = p;
                qq = q;
                ll = lambda;
                psps = psi;
                chichi = chi;
                if aa ~= 0
                    rhoAA = obj.rhoTildFun(aa,q);
                else
                    rhoAA = 0;
                end
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
            med = GenHyperDist.distMedFun(obj.p, obj.lambda, obj.psi, obj.chi, obj.isComplex);
        end
    end
end
