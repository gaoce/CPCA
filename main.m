function main()

% Get data
[data, geneNames] = importData;

% Check data validity
if isempty(data)
    return
end

% Get parameters
[delta, corr] = getParas;

% Check parameter validity
if delta == -1
    return
end

rank = rankCpca(data, geneNames, delta, corr);
exportData(rank);

end