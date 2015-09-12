function [data, geneNames] = importData
% Import data
% Returns:
%   data:  3D data matrix, size(data) = [nTime,nGene,nSampleple]

% Get gene names and gene expression data
[fileNames, pathName, ~] = uigetfile([pwd, '/.gct'], ...
    'Select GCT Files', 'MultiSelect', 'on');

geneNames = {};
for i = 1:length(fileNames)
    path = [pathName, '/', fileNames{i}];
    [names, data(:, :, i)] = gctReader(path); %#ok<AGROW>
    
    % Test if gene names are consistent
    if ~isempty(geneNames)
        if ~all(strcmp(geneNames, names))
            errordlg('Gene Names are not consistent');
            data = [];
            geneNames = {};
            gene2Path = struct();
            % Terminate the function
            return
        else
            continue
        end
    else
        geneNames = names;
    end
end

% % Get Gene Pathway info
% 
% [fileName, pathName, ~] = uigetfile([pwd, '/.gmt'], ...
%     'Select GMT Files');
% 
% gene2Path = gmtReader([pathName, fileName]);

end


function [geneNames, exp] = gctReader(path)
% Read .gct file, get gene expression matrix (gene x time)
% Returns
%   mat: matrix,(time x gene)

fd = fopen(path);

% Skip version line
fgetl(fd);

% Get gene number and time point number
C = textscan(fd, '%d\t%d');
nGene = C{1};
nTime = C{2};

% gene expression matrix
exp = zeros(nGene, nTime);

% gene name cells
geneNames = cell(nGene, 1);

% Skip the header line of the table
fgetl(fd);

% Count the number of tabs in title line, numEx = numTab - 1 + 1
fmtStr = ['%s\t%s\t', repmat('%f\t', 1, nTime-1), '%f\n'];

for i = 1:nGene
    C = textscan(fd, fmtStr);
    
    % Assign gene name
    geneNames{i} = C{1}{1};
    
    % Assign expression values
    exp(i, :) = [C{3:end}];
end

% Transpose the matrix
exp = exp';

end

% function gene2Path = gmtReader(path)
% % Read .gmt file, get gene expression matrix (gene x time)
% % Returns
% %   gene2Path: struct of genes name to cells of pathway names
% 
% gene2Path = struct();
% 
% fid = fopen(path);
% 
% tline = fgetl(fid);
% 
% while ischar(tline)
%     cells = strsplit(tline, '\t');
%     path = cells{1};
%     genes = cells(3:end);
% 
%     for i = 1:length(genes)
%         if isfield(gene2Path, genes{i});
%             n = length(gene2Path.(genes{i}))+1;
%             gene2Path.(genes{i}){n} = path;
%         else
%             gene2Path.(genes{i}) = {path};
%         end
%     end
%     
%     % get next line
%     tline = fgetl(fid);
% end
% 
% fclose(fid);
% end
