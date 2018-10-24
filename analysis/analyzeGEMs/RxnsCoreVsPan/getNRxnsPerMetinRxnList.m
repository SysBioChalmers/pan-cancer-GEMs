function getNRxnsPerMetinRxnList(listfile)

% This functions returns a two column tab-delimited file with the n° of 
% times a metabolite occur in HMR2 (first col) or in the list of rxns (second col)

    %Get nrxns per met in topSub and in HMR
    %Load data
    load ../../data/cHMR3765_20140415.mat
    fid = fopen(listfile,'r');
    t = textscan(fid,'%s','Headerlines',1);
    rxnsTop = t{1};
    cModel = addCompartmentalizedMetNames(cModel);

    %Degree of metabolites (n°rxns per met) in HMR
    degreeHMR = zeros(numel(cModel.mets),1);
    for i = 1:numel(cModel.mets)
        degreeHMR(i) = numel(find(cModel.S(i,:)));
    end

    %Degree of metabolites (n°rxns per met) in HMR
    top.S = cModel.S(:,getIndexes(cModel,rxnsTop,'rxns'));
    degreeTop = zeros(numel(cModel.mets),1);
    for i = 1:numel(cModel.mets)
        degreeTop(i) = numel(find(top.S(i,:)));
    end

    %WrtieMatrix
    writeMatrix([degreeHMR,degreeTop],{'NRxnMetInHMR2';'NRxnMetInList'},cModel.metNamesC,strcat('nrxnmetin_',listfile(1:end-4),'vsHMR2.txt'))
end
function writeMatrix(matrix,header,rownames,filename)
    if nargin < 4
        fid = fopen('mymatrix.txt','wt');
    else
        fid = fopen(filename,'wt');
    end
    if nargin < 3
        rownames = 0;
    end
    if nargin < 2
        header = 0;
    end
    [m,n] = size(matrix);
    if iscell(rownames)
        if iscell(header)
           if length(rownames) ~= m;
            throw(MException('','Rownames must have consistent dimension with matrix'))
           end
           if length(header) ~= n;
            throw(MException('','Header must have consistent dimension with matrix'))
           end
           fprintf(fid,'\t');
           for j = 1:length(header)
               fprintf(fid,'%s\t', header{j});
           end
           fprintf(fid,'\n');
        end
        if isnumeric(matrix) || islogical(matrix)
            for i = 1:size(matrix,1)
                fprintf(fid,'%s\t',rownames{i});
                fprintf(fid,'%4.4f\t',matrix(i,:));
                fprintf(fid,'\n');
            end
            fclose(fid);
        elseif iscell(matrix)
            for i = 1:size(matrix,1)
                fprintf(fid,'%s\t',rownames{i});
                for j = 1:size(matrix,2)
                    fprintf(fid,'%s\t',matrix{i,j});
                end
                fprintf(fid,'\n');
            end
            fclose(fid);
        else
            throw(MException('','Error: input matrix must be either numerical, logical, or cell array of strings'))
        end
    else
        if iscell(header)
           for j = 1:length(header)
               fprintf(fid, '%s\t', header{j});
           end
           fprintf(fid,'\n');
        end
        if isnumeric(matrix)
            for i = 1:size(matrix,1)
                fprintf(fid, '%4.4f\t', matrix(i,:));      
                fprintf(fid,'\n');
            end
            fclose(fid);            
        elseif iscell(matrix)
            for i = 1:size(matrix,1)
                for j = 1:size(matrix,2)
                    fprintf(fid,'%s\t',matrix{i,j});
                end
                fprintf(fid,'\n');
            end
            fclose(fid);
        else
            throw(MException('','Error: input matrix must be either numerical or cell array of strings'))
        end
    end
end
function outModel   = addCompartmentalizedMetNames(inModel)
% addCompartmentalizedMetNames
% This function adds to the model structure a field called metNamesC
% identical to metNames but with the compartment appended at the end.
%
% INPUT
% inModel:  a genome-scale metabolic model structure. Note: if a field
%           called metNamesC is already in the model, it will be
%           overwritten.
% 
% OUTPUT
% outModel:  a genome-scale metabolic model structure with a metNamesC
%            field.
%
% Usage:  outModel   = addCompartmentalizedMetNames(inModel)
% 2013-07-17 Francesco Gatto
    outModel = inModel;
    if isfield(outModel,'metNamesC')
        outModel = rmfield(outModel,'metNamesC');
    end
    mets                = outModel.metNames;
    nMets               = length(mets);
    comps               = outModel.comps(outModel.metComps);
    metNamesC           = cellfun(@(a,b,c,d) [a,b,c,d],mets,repmat({'['},nMets,1),comps,repmat({']'},nMets,1),'uni',false);
    outModel.metNamesC  = metNamesC;
end