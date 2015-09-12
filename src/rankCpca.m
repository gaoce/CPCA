function [rank] = rankCpca(data, geneNames, delta, corr)
% RANKCPCA implement CPCA algorithm, which rank genes based on their
% temporal variations across different samples
% Argument:
%     data:  3D data matrix, size(data) = [nTime,nGene,nSampleple]
%     feat:  contain gene names (feat.gene) and correponding pathway 
%            (feat.pathway)
%     delta: percentage of information to retain in PCA
%     corr:  logical, 1 correlation matrix (default), 0 covariance

% Author:   Ce Gao
% Created:  08-20-2012
% Revised:  08-31-2012
%           11-15-2012 => refactoring and documentation
%           12-10-2012 => documentation
%           03-27-2013 => documentation
%           09-12-2015 => Modified functionality, and added to CPCA toolbox
% Toolbox:  CPCA, based on MATLAB R2015a
% Reference:
% 1. Between-groups comparison of principal components
% 2. Feature subset selection and feature ranking for multivariate time
%    series


DCPC = computeDCPC(data, delta, corr);
rank = rankGene(DCPC, geneNames);

end % primary function

function DCPC = computeDCPC(data, delta, corr)
% COMOUTEDCPC calculates  descriptive common principal components (DCPCs)
% of a multidimensional matrix data, note the dimension order is strict
% Argument:
%     data:  3D data matrix, size(data) = [nTime,nGene,nSampleple]
%     delta: percentage of information to retain in PCA
%     corr:  logical, 1 correlation matrix (default), 0 covariance
% Output:
%     DCPC:  size(DCPC) = [p,nGene] (p depends on delta). DCPC(i,j) is the
%            common score on ith princomp for jth gene

%% Argument passing
if nargin < 3
    corr = 1;
end
[~,nGene,nSample] = size(data);

%% PCA to each items using SVD
loading = zeros(nGene,nGene,nSample);
nVar    = zeros(nSample,1);

switch corr
    case 1 % zscore
        dataNorm = zscore(data);
    case 0 %just center
        dataNorm = bsxfun(@minus,data,mean(data));
    otherwise
        error('Incorrect argument!');
end
dataNorm(isnan(dataNorm)) = 0;

for i = 1:nSample
    [~,S,V]        = svd(dataNorm(:,:,i));
    loading(:,:,i) = V';  %in loadings, cols are real genes
    percentVar     = 100*cumsum(diag(S))/sum(diag(S));
    nVar(i)        = sum(percentVar < delta) + 1;
end
p = max(nVar);
%% Compute common PCA
% Note: in loadings(:,:,i), each row containing coefficients for one
% principal component. The rows are in order of decreasing component
% variance

H = zeros(nGene,nGene);
L = loading(1:p,:,:);
for i = 1:nSample
    H = H + L(:,:,i)'*L(:,:,i);
end

[~,~,VH] = svd(H);
DCPC = VH(:,1:p)';
end

function rank = rankGene(DCPC, geneNames)
% rankGene rank genes according their DCPC scores
% Argument:
%     DCPC:  size(DCPC) = [p,nGene] (p depends on delta). DCPC(i,j) is the
%            common score on ith princomp for jth gene
% Output:
%     rank.gene:  gene names, with decreasing scores
%     rank.path:  gene pathway, with decresing scores
%     rank.score: scores

score = sum(DCPC.^2,1)';
[~,IDX] = sort(score,'descend');

%% output wrapper
rank.gene  = geneNames(IDX);
rank.score = score(IDX);
end

