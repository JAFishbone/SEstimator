% Copyright (C) 2020-2022 Justin A Fishbone
classdef SEstimator < handle % S-Estimator

    % Required constant subclass properties
    properties (Abstract, SetAccess = protected, GetAccess = public)
        isSq(1,1)        logical         % Is this an S-q estimator?
        maxTun(1,1)      double          % Maximum tuning parameter value to use for enlargement logic above 
    end                                  % (must result in: numerically computable result, proper S-estimator)
        
    % Public properties
    properties
        tun0(1,1)        double          % Initial tuning parameter
        tun(1,1)         double          % Current tuning parameter
        p(1,1)           double = NaN;   % Data dimension
        pEff(1,1)        double = NaN;   % Distribution's effective p (i.e. 2p for complex-valued case)
        isComplex(1,1)   logical         % Is this a complex-valued case?
        b                double          % S-estimator parameter b; didn't force size (1,1) so property can be empty
        minPropVld(1,1)  double = 1.5;   % Proportion of p equal to the minimum number of samples with positive weights
        tol              double = 1e-6;  % Solution convergence criteria
    end
        
    % Private properties
    properties (Access = protected)
        tunEnlarged(1,1) double = NaN;   % Enlarged tuning parameter (if needed)
    end
    
    % Subclass methods
    methods (Abstract)
        rho = rhoTildFun(obj,d);            % Decreasing rho (rho-tilde) function
        w = wTildFun(obj,d);                % w-tilde function
        [a,c,rhoA,rhoC] = extremesFun(obj); % Estimator extremes function
        med = distMedFun(obj);              % Returns median of expected PDF
    end
    
    % Protected methods
    methods (Access = protected)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% Constructor
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function obj = SEstimator(varargin)
        % obj = SEstimator(tun0, tol, b, isComplex) -- All optional inputs
            if nargin > 0
                obj.tun0        = varargin{1};
                obj.tun         = varargin{1};
            end
            if nargin > 1 && ~isempty(varargin{2})
                obj.tol         = varargin{2};
            end
            if nargin > 2
                obj.b           = varargin{3};
            end
            if nargin > 3
                obj.isComplex   = varargin{4};
            end
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% estimateScale(): Determine sigHat -- Used by obj.estimate()
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function sigHat = estimateScale(obj, d)
            % Uses method discussed by Maronna et al. 2019 Section 9.4.2, Pg 368
            % Find lower bound of search
            s0 = median(d);
            lb = s0;
            while (mean(obj.rhoFun(d/lb)) - obj.b) < 0
                lb = lb/2;
            end
            % Find upper bound of search
            ub = s0;
            while (mean(obj.rhoFun(d/ub)) - obj.b) > 0 && ~isinf(ub)
                ub = ub*2;
            end
            % Solve the S-estimator equation for sigHat
            x0 = [lb, ub];
            sigHat = fzero( (@(s)mean(obj.rhoFun(d./s)) - obj.b), x0 );
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% numValidSamples(): Number of samples with positive weights -- Used by obj.estimate()
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function numValid = numValidSamples(obj, ds, tun)
            if nargin < 3
                tun = obj.tun;
            end
            [a,c,~,~] = obj.extremesFun(tun);
            numValid = sum(ds > a & ds < c);
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% enlargeTun(): Enlarge the tuning parameter so sufficient samples have positive weight -- Used by obj.estimate()
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function enlargeTun(obj, ds, minValidSamples)
            if obj.numValidSamples(ds, obj.maxTun) < minValidSamples
                obj.tunEnlarged = obj.maxTun;
            else
                obj.tunEnlarged = fzero(@(tuning) obj.numValidSamples(ds, tuning) - minValidSamples, [obj.tun0, obj.maxTun]);
            end
            obj.tun = obj.tunEnlarged;
        end
    end
        
    % Public methods
    methods
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% rhoFun(): Rho function
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function rho = rhoFun(obj, t)
            % Decreasing rho (rho-tilde) function
            if obj.isSq && obj.isComplex 
            % Complex-valued S-q case
                rho = 2^(1-obj.tun) * obj.rhoTildFun(2*t);
            else
                rho = obj.rhoTildFun(t);
            end
            
            % Normalize (unless S-q with q=1)
            if ~(obj.isSq && obj.tun == 1)
                % Get the extreme points
                [a,c,rhoA,rhoC] = obj.extremesFun;
                if obj.isSq && obj.isComplex
                % Correction for complex-valued case
                   a = a/2;
                   c = c/2;
                   rhoA = rhoA*2^(1-obj.tun);
                   rhoC = rhoC*2^(1-obj.tun);
                end
                % Shift and scale rho
                rho = (rho-rhoA)/(rhoC-rhoA);
                % Make non-decreasing
                rho(t<=a) = 0;
                rho(t>=c) = 1;
            end
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% wFun(): Weight function
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function w = wFun(obj, t)
            % Weight function corresponding to decreasing rho-function (rho-tilde)
            if obj.isSq && obj.isComplex 
            % Complex-valued S-q case
                w = 2^(2-obj.tun) * obj.wTildFun(2*t);
            else
                w = obj.wTildFun(t);
            end
            
            % Set zero weights unless S-q with q=1
            if ~(obj.isSq && obj.tun == 1)
                % Get the extreme points
                [a,c,~,~] = obj.extremesFun;
                if obj.isSq && obj.isComplex
                % Correction for complex-valued case
                   a = a/2;
                   c = c/2;
                end
                % Set weights to 0 outside points a & c
                w(t<=a) = 0;
                w(t>=c) = 0;
            end
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% estimate(): Main estimation function 
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [locHat, shapeHat, scatterHat, sigHat, status] = estimate(obj, X, loc0, shape0, tun, tol, b)
        % [locHat, shapeHat, scatterHat, sigHat, status] = estimate(obj, X, loc0, shape0)
        % [locHat, shapeHat, scatterHat, sigHat, status] = estimate(obj, X, loc0, shape0, tun, tol, b)
            
            %--------------------------------------------------------------
            %%% PROCESS INPUTS
            %--------------------------------------------------------------
            
            % Determine sample size
            n = size(X,2);
            if obj.p ~= size(X,1) % Reset obj.p if necessary
                obj.p = size(X,1);
            end
            
            % Update whether this is a complex-vlaued case
            if isreal(X)
                obj.isComplex = false;
                obj.pEff      = obj.p;
            else
                obj.isComplex = true;
                obj.pEff      = 2*obj.p;
            end
            
            % Tuning parameter
            if nargin > 4 && ~isempty(tun)
                obj.tun0 = tun;
                obj.tun  = tun;
            end
            
            % Numerical convergence tolerance/threshold
            if nargin > 5 && ~isempty(tol)
                obj.tol  = tol;
            end
            
            % Breakdown point parameter
            if nargin > 6 % Will assign empty if input is empty so we can recalculate be below
                obj.b    = b;
            end
            if isempty(obj.b)
                obj.b = 0.5-0.5*obj.p/n; % Lopuhaa 1991
            end

            % Initial estimates
            shapeHatPrev = shape0/det(shape0)^(1/obj.p); % Ensure true shape matrix
            locHatPrev   = loc0;

            %--------------------------------------------------------------
            %%% PERFORM ESTIMATE
            %--------------------------------------------------------------
            i = 1;
            % Calculate square Mahalanobis distance and scale estimate
            d = mdist(X,locHatPrev,shapeHatPrev);
            sigHat(i) = obj.estimateScale(d);
            % Find tun to ensure that at least p*minPropValid elements have w>0
            ds=d/sigHat(i);
            if obj.numValidSamples(ds) < (ceil(obj.p*obj.minPropVld) + 1) && obj.tun0 ~= 1
                % Enlarging tuning parameter for initial iteration because too many samples have 0 weight
                enlarged_tun = true;
                obj.enlargeTun(ds, ceil(obj.p*obj.minPropVld) + 1);
                sigHatPrev = obj.estimateScale(d);
            else
                enlarged_tun = false;
                sigHatPrev = sigHat(i);
            end
            status.enlarged_initial_tun = enlarged_tun;

            holdon = false;
            err = inf;
            while err > obj.tol
                % Weights
                w = obj.wFun(d./sigHat(i));
                % See if tun needs to be enlarged for 2nd iteration, if so, hold it for all iterations
                if enlarged_tun && sum(w~=0) < ceil(obj.minPropVld*obj.p) && ~holdon
                    % Holding tuning parameter
                    obj.tun = obj.tunEnlarged;
                    holdon = true;
                    w = obj.wFun(d./sigHat(i));
                    sigHatPrev = obj.estimateScale(d);
                end
                if sum(w~=0) < 1+obj.p
                % This condition is very rare, but will result in a singular shape matrix
                % To Do: Update algorithm to handle this case by inflating tuning parameter
                    error('Reached point where there is an insufficienct number of positive-weight samples');
                end

                % Advance to next iteration
                i = i + 1; 

                % Estimate location
                locHatNext = sum(w.*X, 2) ./ sum(w);

                % Estimate shape matrix
                for j = 1:n
                    tmpS0(:,:,j) = (X(:,j) - locHatNext) * (X(:,j) - locHatNext)'; %#ok<AGROW>
                end
                tmpW = repmat( reshape(w, 1,1,n), [obj.p,obj.p,1] );
                tmpS = mean( tmpW .* tmpS0, 3);
                while det(tmpS) > 1e100 % Sometimes c can get too large and det(c) computes to inf in MATLAB, so shrink it
                    tmpS = tmpS/10;
                end
                shapeHatNext = tmpS ./ det(tmpS)^(1/obj.p);

                % Calculate Mahalanobis distance and scale estimate
                d = mdist(X,locHatNext,shapeHatNext);
                sigHat(i) = obj.estimateScale(d);

                % Can sometimes overstep, which increases sigHat
                % If so, take partial step using method from Maronna et al. 2019 Section 9.6.3, Pg 371
                xi0 = 0.7;
                xi = xi0;
                while sigHat(i) >= sigHatPrev && xi > 0.001 % Eq 9.25, Pg 371 

                    % Updated location and shape estimates
                    locHatNext = (1-xi)*locHatPrev + xi*locHatNext; % Eq 9.26, Pg 371
                    shapeHatNext = (1-xi)*shapeHatPrev + xi*shapeHatNext; % Eq 9.26, Pg 371
                    shapeHatNext = shapeHatNext ./ det(shapeHatNext)^(1/obj.p);

                    % Calculate Mahalanobis distance and scale estimate
                    d = mdist(X,locHatNext,shapeHatNext);
                    sigHat(i) = obj.estimateScale(d);

                    % Shrink xi in case another iteration is needed
                    xi = xi0*xi;
                    if xi < 0.001
                        status.hit_xi_limit = true;
                    end
                end

                % Calculate error metric
                % If we've enlarged tun, then ensure we've had 3 iterations
                if ~enlarged_tun || i > 3 % Set to 3 b/c if we hit xi limit on i==3
                   err = trace(shapeHatNext\shapeHatPrev)-obj.p; 
                end

                % Reset tun if we're not holding it for the run
                if ~holdon && i==2
                    obj.tun = obj.tun0;
                    sigHatPrev = obj.estimateScale(d);
                else
                    sigHatPrev = sigHat(i);
                end

                % Prepare for next iteration
                locHatPrev   = locHatNext;
                shapeHatPrev = shapeHatNext;
            end

            % Grab shape->scatter scaling factor
            % These methods fairly (but not the most) robust, so catch errors and assign NaN
            warning('off');
            try
                med = obj.distMedFun();
            catch
                med = NaN;
            end
            warning('on');
            
            % Set output variables
            status.hold_tun = holdon;
            status.final_tun = obj.tun;
            locHat = locHatNext;
            shapeHat = shapeHatNext;
            % Scale the scatter matrix
            if ~isnan(med)
                % Scale the scatter matrix by using the median of the shape-based squared Mahalanobis distances
                scatterHat = shapeHat*median(d)/med; % Estimate the scatter matrix
                status.scatterScaleType = 'MEDIAN';
            else
                % Scale the scatter matrix by using sigma-hat and the S-estimator rho function constraint
                scatterHat = shapeHat*sigHat(end);
                status.scatterScaleType = 'RHO';
            end
        end
    end
end