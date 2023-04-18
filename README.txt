Copyright (C) 2020-2023 Justin A Fishbone

This package contains the S-Estimator package that includes the Sq-estimator discussed in:

Fishbone, J., Mili, L. (Accepted 2022). "New Highly Efficient High-Breakdown Estimator of Multivariate Scatter and Location for Elliptical Distributions." Canadian Journal of Statistics, doi: 10.1002/cjs.11770
https://onlinelibrary.wiley.com/doi/10.1002/cjs.11770

Fishbone, J.A., Mili, L. (2023). "Highly Robust Complex Covariance Estimators with Applications to Sensor Array Processing." IEEE Open Journal of Signal Processing, doi: 10.1109/OJSP.2023.3261806
https://ieeexplore.ieee.org/document/10081068


INSTALLATION
1) Run add_paths from this directory


EXECUTION
Execute an S-estimate using the sEst command:
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
% 'QCAUCHY'       - Sq Cauchy                   - q     - N/A
% 'QGAUSS'        - Sq Gaussian (i.e. Normal)   - q     - N/A
% 'QGAUSSIAN'     - SAME AS 'QGAUSS'
% 'QGENHYPER'     - Sq Generalized Hyperbolic   - q     - distParms.chi, distParms.lambda, distParms.psi
% 'QKOTZ'         - Sq Kotz Type                - q     - distParms.N, distParms.r, distParms.s
% 'QMVHYPERBOLIC' - Sq Multivariate Hyperbolic  - q     - distParms.chi, distParms.psi
% 'QNORMAL'       - SAME AS 'QGAUSS'
% 'QNORMINVGAUSS' - Sq Normal Inverse Gaussian  - q     - distParms.chi, distParms.psi
% 'QLAPLACE'      - Sq Laplace                  - q     - N/A
% 'QLOGISTIC'     - Sq Logistic                 - q     - N/A
% 'QPEARSONII'    - Sq Pearson Type II          - q     - distParms.m
% 'QPEARSONVII'   - Sq Pearson Type VII         - q     - distParms.N, distParms.s
% 'QT'            - Sq t-Distribution           - q     - distParms.nu
% 'QVARGAMMA'     - Sq Variance Gamma           - q     - distParms.lambda, distParms.psi
%


TO ADD ADDITIONAL S-ESTIMATORS TO sEst SUB-PACKAGE
1a) Copy and rename one of the existing S-Estimator subclasses in
    .\sEst\SEstimator\subclasses\
Recommend using either sRocke.m or qKotz.m.

1b) Update all properties and methods for the new estimator

2) Update subfunction initEst() in sEst.m with the new case for the new estimator


TO ADD ADDITIONAL DISTRIBUTION CLASSES TO Dist SUB-PACKAGE
1) Copy and rename one of the existing distribution classes in
    .\Dist\
Recommend using GaussDist.m for a simple example.
Recommend using KotzDist.m for an example without closed-form expressions.