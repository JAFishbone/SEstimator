% Copyright (C) 2020-2022 Justin A Fishbone
function X = randEllip(f, p, n, mu, sigma, R_REF_MAX_SURV, ABS_TOL, REL_TOL, TOL, N_REF_SAMPLES, IS_COMPLEX)
% X = randEllip(f, p, n, mu, sigma)
% X = randEllip(..., R_REF_MAX_SURV, ABS_TOL, REL_TOL, TOL, N_REF_SAMPLES)
% X = randEllip(..., R_REF_MAX_SURV, ABS_TOL, REL_TOL, TOL, N_REF_SAMPLES, IS_COMPLEX)
% Generate multivariate elliptically symmetric random vectors X=[x1, ..., xn]
%
% Method uses random direction times a random magnitude
% Random magnitude is determined using a sampled CDF function stored in memory
%
% INPUTS:
% f              - Function handle for the PDF of the absolute Mahalanobis distance (not squared Mahalanobis distance)
% p              - Dimension
% n              - Number of random vectors
% mu             - (Optional) px1 location vector (default is zeros)
% sigma          - (Optional) pxp scatter matrix (default is identity)
% R_REF_MAX_SURV - (Optional) Vector of point to calculate CDF at to find the approximate max (1-TOL) of the CDF
% ABS_TOL        - (Optional) Absolute tolerance used for integral function
% REL_TOL        - (Optional) Relative tolerance used for integral function
% TOL            - (Optional) Tolerance that defines the numerical start (i.e. not 0) and end (i.e. not 1) of the CDF
%                  (this sets the numerical precision at the tails of the CDF)
% N_REF_SAMPLES  - (Optional) Number of sample points to calculate the CDF. 
%                  These samples are used to interpolate the final random magnitude
% IS_COMPLEX     - (Optional) Whether to generate complex-valued data (default is false)
%
% OUTPUTS:
% X - pxn matrix of n p-dimensional elliptically distributed random vectors
%
% NOTES:
% Function handle f cannot have nested function handles!
%

% Process inputs/set defaults
if nargin < 6 || isempty(R_REF_MAX_SURV)
    R_REF_MAX_SURV = 2.^(0:45);
end
if nargin < 7 || isempty(ABS_TOL)
    ABS_TOL = 1e-12;
end
if nargin < 8 || isempty(REL_TOL)
    REL_TOL = 1e-12;
end
if nargin < 9 || isempty(TOL)
    TOL = 1e-10;
end
if nargin < 10 || isempty(N_REF_SAMPLES)
    N_REF_SAMPLES = 1e4;
end
if nargin < 11 || isempty(IS_COMPLEX)
    IS_COMPLEX = false;
end

% Grab details about the function handle
f_details = functions(f); % This could break in a future release of MATLAB

% Since CDF is not known in closed form, calculate a reference CDF with which to interpolate
% User persistent variables so sampled reference function only has to be calculated once
persistent r_ref FRef f_func f_ws pp pR_REF_MAX_SURV pABS_TOL pREL_TOL pTOL pN_REF_SAMPLES;
if isempty(r_ref) || pp ~= p || ~strcmp(f_details.function, f_func) || ~isequal(f_details.workspace{1}, f_ws) ...
                  || ~isequal(pR_REF_MAX_SURV, R_REF_MAX_SURV) || pABS_TOL ~= ABS_TOL || pREL_TOL ~= REL_TOL || pTOL ~= TOL || pN_REF_SAMPLES ~= N_REF_SAMPLES % Can probably get rid of these checks
    % Temporarily disable this warning
    warning('off', 'MATLAB:integral:NonFiniteValue');

    % Integrate PDF at sampled survey values to find approximate max (1-TOL) of CDF
    f_ref_max_surv = arrayfun(@(r) integral(@(rr) f(rr), 0, r, 'reltol', REL_TOL, 'abstol', ABS_TOL), R_REF_MAX_SURV);
    f_ref_max_surv(isinf(f_ref_max_surv)) = NaN; % Inf is a numerical error
    r_ref_max = find(f_ref_max_surv > 1-TOL, 1, 'first');
    if isempty(r_ref_max)
        error('Unable to find max of CDF');
    end
    
    % Reference Mahalanobis distance for CDF F(r)
    r_ref = linspace(0, R_REF_MAX_SURV(r_ref_max), N_REF_SAMPLES);
    % Integrate PDF to get CDF at each reference value r_ref
    FRef = arrayfun(@(r) integral(@(rr) f(rr), 0, r, 'reltol', REL_TOL, 'abstol', ABS_TOL), r_ref(2:end));
    % Inf is a numerical error
    FRef(isinf(FRef)) = NaN;
    % Add F(0)=0
    FRef = [0, FRef];
    
    % Trim the ends to remove any noise
    lastZeroInd = find(FRef<TOL, 1, 'last');
    if isempty(lastZeroInd)
        lastZeroInd = 1;
    end
    firstOneInd = find(FRef>1-TOL, 1, 'first');
    if isempty(firstOneInd)
        firstOneInd = length(FRef);
    end
    r_ref = r_ref(lastZeroInd:firstOneInd);
    FRef  =  FRef(lastZeroInd:firstOneInd);
    
    % Ensure there's a zero point
    if FRef(1) ~= 0
        FRef  = [0, FRef];
        r_ref = [0, r_ref];
    end
    
    % Save persistent variables
    f_func = f_details.function;
    f_ws   = f_details.workspace{1};
    pp     = p;
    pR_REF_MAX_SURV = R_REF_MAX_SURV;
    pABS_TOL = ABS_TOL;
    pREL_TOL = REL_TOL;
    pTOL = TOL;
    pN_REF_SAMPLES = N_REF_SAMPLES;
    
    % Re-enable this warning
    warning('on', 'MATLAB:integral:NonFiniteValue');
end

% Determine random magnitude
% Uniform random values
u = rand(1,n);
% Random magnitudes
r = interp1( FRef, r_ref, u, 'pchip', NaN );
if any(isnan(r))
    warning('Random value is above support');
    r(isnan(r)) = r_ref(end);
end
% Warn if below support
if any(u < TOL)
    warning('Random value is below support');
end

% Determine random direction
% Generate random vectors on the unit p-sphere
if IS_COMPLEX
    uSphTmp = complex(randn(p,n), randn(p,n)); % Don't need /sqrt(2) b/c we normalize below
else
    uSphTmp = randn(p,n);
end
uSph = uSphTmp ./ vecnorm(uSphTmp,2,1);

% Standard multivariate logistic random variable
X = r.*uSph;

% Multiply scatter
if nargin < 5 || isempty(sigma) || isequal(sigma, eye(p))
    s = 1;
else
    s = cholcov(sigma)';
end
X = s*X;

% Add mean
if nargin > 3 && ~isempty(mu)
    X = mu + X;
end
end
