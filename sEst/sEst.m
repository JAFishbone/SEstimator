% Copyright (C) 2020-2022 Justin A Fishbone
function [locHat, shapeHat, scatterHat, sigHat, status] = sEst(X, estFun, distParms, tun0, loc0, shape0, varargin)
% [locHat, shapeHat, scatterHat, sigHat, status] = sEst(X, estimator, distParms, tun, loc0, shape0)
% [locHat, shapeHat, scatterHat, sigHat, status] = sEst(..., tol, b)
% S-estimator of location, shape, and scatter matrices
%
% INPUTS:
% X          - pxn matrix of observed data samples
% estFun     - String specifying which S-estimator - See available options below
% distParms  - Structure of additional estimator parameters related to underlying distribution - See details of estimators below
%              (when N/A, input empty value [])
% tun        - Estimator tuning parameter (not applicable to S-Bisq)
% loc0       - Initial guess for location vector (px1 vector)
% shape0     - Initial guess of shape (or scatter) matrix (pxp matrix)
% tol        - (Optional) Convergence criterion (scalar)
% b          - (optional) S-estimator b value in equation mean(rho(.))=b (scalar) (default is maximum breakdown value)
%
% OUTPUTS:
% locHat     - Location estimate (px1 vector)
% shapeHat   - Shape matrix estimate (pxp matrix)
% scatterHat - Scatter matrix estimate (pxp matrix) (For S-Bisq & S-Rocke, this is only appropriate if underlying data is Gaussian)
% sigHat     - Vector of sigma estimates for each iteration
% status     - Status information from the estimation
%
% DETAILS ON ESTIMATOR FUNCTIONS
% estFun          FUNCTION NAME                 tun     REQUIRED distParms
% 'BISQ'          - S-Bisquare                  - N/A   - N/A
% 'ROCKE'         - S-Rocke                     - gamma - N/A
% 'QCAUCHY'       - S-q Cauchy                  - q     - N/A
% 'QGAUSS'        - S-q Gaussian (i.e. Normal)  - q     - N/A
% 'QGAUSSIAN'     - SAME AS 'QGAUSS'
% 'QGENHYPER'     - S-q Generalized Hyperbolic  - q     - distParms.chi, distParms.lambda, distParms.psi
% 'QK'            - S-q K-distribution          - q     - distParms.nu
% 'QKOTZ'         - S-q Kotz Type               - q     - distParms.N, distParms.r, distParms.s
% 'QLAPLACE'      - S-q Laplace                 - q     - N/A
% 'QLOGISTIC'     - S-q Logistic                - q     - N/A
% 'QMVHYPERBOLIC' - S-q Multivariate Hyperbolic - q     - distParms.chi, distParms.psi
% 'QNORMAL'       - SAME AS 'QGAUSS'
% 'QNORMINVGAUSS' - S-q Normal Inverse Gaussian - q     - distParms.chi, distParms.psi
% 'QPEARSONII'    - S-q Pearson Type II         - q     - distParms.m
% 'QPEARSONVII'   - S-q Pearson Type VII        - q     - distParms.N, distParms.s
% 'QT'            - S-q t-Distribution          - q     - distParms.nu
% 'QVARGAMMA'     - S-q Variance Gamma          - q     - distParms.lambda, distParms.psi
% 'QW       '     - S-q W-distribution          - q     - distParms.s, distParms.b
%

% Created by XX_BLINDEDNAME_XX 2021
% Developed and tested in MATLAB R2019a

% Initialize estimator object
p = size(X,1);
est = initEst(estFun, p, distParms);

% Calculate estimate
[locHat, shapeHat, scatterHat, sigHat, status] = est.estimate(X, loc0, shape0, tun0, varargin{:});
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% SUPPORT FUNCTIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function est = initEst(estFun, p, distParms)
% Initialize the estimator functions
%

% Identify the estimator functions. Try to catch common alternate specifiers
switch upper(estFun)
    case {'BISQ', 'BISQUARE'}
        est = sBisq();
    case {'ROCKE', 'ROCK'}
        est = sRocke();
    case 'QCAUCHY'
        est = qPearsonVII(1, (1+p)/2);
    case {'QGAUSS', 'QGAUSSIAN', 'QNORMAL'}
        est = qGauss();
    case {'QGENHYPER', 'QGENERALIZEDHYPERBOLIC'}
        est = qGenHyper(distParms.lambda, distParms.psi, distParms.chi);
    case {'QK', 'QKDIST', 'QKDISTRIBUTION'}
        est = qGenHyper(distParms.nu, 2*distParms.nu, 0);
    case {'QKOTZ', 'QKOTZTYPE'}
        est = qKotz(distParms.s,distParms.N,distParms.r);
    case {'QMVHYPER', 'QMVHYPERBOLIC', 'QMULTIVARHYPER', 'QMULTIVARHYPERBOLIC'}
        est = qGenHyper((p+1)/2,          distParms.psi, distParms.chi);
    case 'QNORMINVGAUSS'
        est = qGenHyper(   -1/2,          distParms.psi, distParms.chi);
    case 'QLAPLACE'
        est = qGenHyper(      1,                       2,            0);
    case 'QLOGISTIC'
        est = qLogistic();
    case {'PEASRONII', 'PEARSONTYPEII', 'PEASRON2', 'PEARSONTYPE2'}
        est = qPearsonII(distParms.m);
    case {'PEASRONVII', 'PEARSONTYPEVII', 'PEASRON7', 'PEARSONTYPE7'}
        est = qPearsonVII(distParms.s, distParms.N);
    case {'QT', 'QSTUDENT', 'QSTUDENTT'}
        est = qPearsonVII(distParms.nu, (distParms.nu + p)/2);
    case {'QVARGAMMA', 'QVARIANCEGAMMA'}
        est = qGenHyper(distParms.lambda, distParms.psi,             0);
    case {'QW', 'QWDIST', 'QWDISTRIBUTION'}
        est = qKotz(distParms.s, distParms.s-1, 1/(2^distParms.s * distParms.b));
    otherwise
        fprintf(2, 'ERROR: The value ''%s'' is an unexpected value for estFun\n', estFun);
        fprintf(2, 'Allowed values are:\n\n');
        fprintf(2, 'estFun          FUNCTION NAME                 tun     REQUIRED distParms\n');
        fprintf(2, '''BISQ''          - S-Bisquare                  - N/A   - N/A\n');
        fprintf(2, '''ROCKE''         - S-Rocke                     - gamma - distParms.alpha (an alternative to gamma, overrides gamma)\n');
        fprintf(2, '''QCAUCHY''       - S-q Cauchy                  - q     - N/A\n');
        fprintf(2, '''QGAUSS''        - S-q Gaussian (i.e. Normal)  - q     - N/A\n');
        fprintf(2, '''QGENHYPER''     - S-q Generalized Hyperbolic  - q     - distParms.chi, distParms.lambda, distParms.psi\n');
        fprintf(2, '''QK''            - S-q K-Distribution          - q     - distParms.nu\n');
        fprintf(2, '''QKOTZ''         - S-q Kotz Type               - q     - distParms.N, distParms.r, distParms.s\n');
        fprintf(2, '''QMVHYPERBOLIC'' - S-q Multivariate Hyperbolic - q     - distParms.chi, distParms.psi\n');
        fprintf(2, '''QNORMINVGAUSS'' - S-q Normal Inverse Gaussian - q     - distParms.chi, distParms.psi\n');
        fprintf(2, '''QLAPLACE''      - S-q Laplace                 - q     - N/A\n');
        fprintf(2, '''QLOGISTIC''     - S-q Logistic                - q     - N/A\n');
        fprintf(2, '''QPEARSONII''    - S-q Pearson Type II         - q     - distParms.m\n');
        fprintf(2, '''QPEARSONVII''   - S-q Pearson Type VII        - q     - distParms.N, distParms.s\n');
        fprintf(2, '''QT''            - S-q t-Distribution          - q     - distParms.nu\n');
        fprintf(2, '''QVARGAMMA''     - S-q Variance Gamma          - q     - distParms.lambda, distParms.psi\n');
        fprintf(2, '''QW''            - S-q W-Distribution          - q     - distParms.s, distParms.b\n\n');
        error('Unexpected value for estFun.');
end
end