% Copyright (C) 2020-2022 Justin A Fishbone
function [muHat, omegaHat, sigmaHat] = mve(X, nSubSampling)
% function [muHat, omegaHat, sigmaHat] = mve(X, nSubSampling)
% Minimum volume ellipsoid location and scatter estimates
% following Maronna et al. 2019 Section 6.8.3
%
% INPUTS:
% X           - n x p matrix of samples
% nSubSamples - (optional) number of subsample sets to use
%
% OUTPUTS:
% muHat       - Location estimate
% omegaHat    - Shape matrix estimate
% sigmaHat    - Scatter matrix estimate

% For reproducibility, seed random number generator
rngState = rng;
rng(0);
% Ignore matrix singularity warnings, this is expected with small subsamples
warning('off', 'MATLAB:nearlySingularMatrix');

% Set defaults
if nargin < 2 || isempty(nSubSampling)
    nSubSampling = 5000;
end

% Grab parameters of interest
[p,n] = size(X);

% Determine subsampling indices
k = p+1;
inds = random_subsampling(n, k, nSubSampling);

% Process each subsample, and keep track of the one with the smallest median distance
smallestMedDist = inf;
for j = 1:nSubSampling
    XSub           = X(:,inds(:,j));
    muHatTmp       = mean(XSub,2);
    XSubCent       = (XSub-muHatTmp);
    sigmaHatTmp    = XSubCent*XSubCent'/k;
    omegaHatTmp       = sigmaHatTmp/det(sigmaHatTmp)^(1/p);
    % Distances of all samples based on the location and shape of the subsample
    dists = mdist(X,muHatTmp,omegaHatTmp);
    % Re-calculate the estimates using the half of the samples with smallest distances
    XEnlarg = X( :, dists <= median(dists) );
    muHatTmp       = mean(XEnlarg,2);
    XSubCent       = (XEnlarg-muHatTmp);
    sigmaHatTmp    = XSubCent*XSubCent'/size(XEnlarg,2);
    omegaHatTmp       = sigmaHatTmp/det(sigmaHatTmp)^(1/p);
    % Distances of all samples based on the updated location and shape of the subsample    
    dists = mdist(X,muHatTmp,omegaHatTmp);
    % Does this subsample have the smallest median distance?
    medDist = median(dists);
    if medDist < smallestMedDist
        % This is out current estimate
        muHat    = muHatTmp;
        sigmaHat = sigmaHatTmp;
        omegaHat = omegaHatTmp;
        smallestMedDist = medDist;
    end
end

% Return random number generator to previous state
rng(rngState);

% Turn warning back on
warning('on', 'MATLAB:nearlySingularMatrix');
end



function inds = random_subsampling(n, k, nSubSampling)
% function myUniShufs = random_subsampling(n, k, nSubSampling)
% Generate nSubSamples subsamples of size k from a total sample size of length n
%
% Adapted from:
% https://www.mathworks.com/matlabcentral/answers/230988-how-to-return-x-number-of-unique-subsets-combinations-of-n-numbers-taken-k-at-a-time
%
% Inputs:
% n - Total number of samples
% k - Subsample size
% nSubSamples - number of subsamples of size k from the n total samples
%

mixer = false(n,1);
mixer(1:k) = true;
inds = zeros(k,nSubSampling);
% Randomly shuffle with margin
for i = 1:nSubSampling*10 % Tune this for efficiency
   mixer = mixer(randperm(n)); 
   inds(:,i) = find(mixer);
end
% Remove duplicate entries
inds = unique(inds', 'stable', 'rows')';
while size(inds,2) < nSubSampling % Recursion here ensures we have enough unique subsamplings
    inds = random_subsampling(n, k, nSubSampling);
end
% Return the requested number of subsamples
inds = inds(:,1:nSubSampling);

end



















