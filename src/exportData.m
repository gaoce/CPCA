function exportData(rank)
% EXPORTDATA writse ranking and clustering result to output files, in TSV 
% (tab separated value) format. 
% 
% Input:
%     fileName: file name to write
%     rank:     ranking result with 3 fields: gene, score, pathway
% Output:
%    for rank:
%        geneName    score    pathway

% Authors:   Ce Gao
% Created:  11-15-2012
% Revised:  04-12-2012 => get rid of clustering part
% Toolbox:  CPCA, based on MATLAB R2012a

[fileName, pathName, ~] = uiputfile([pwd, '/.rnk'], ...
    'Enter Gene List File Name:', 'cpca.rnk');

fh = fopen([pathName, fileName],'w');
for j  = 1:length(rank.gene)
    fprintf(fh,'%s\t%f\n',rank.gene{j},rank.score(j));
end

fclose(fh);