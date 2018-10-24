function [elInclMat,elRef,PSID] = getElementInclusionMatrix(type,status,PSID,verbose)
%
% getElementInclusionMatrix
%
% It returns a binary matrix if an element (gene, rxn, or met) from a reference model is included
% in the models specified by PSID or all (by default)
%
% INPUT
% type : a string that can be "genes", "rxns", or "mets".
% status: a string that can be "T" (def.) or "N", if models to be loaded
%         are from tumor or normal samples.
% PSID : a cell array of strings containing the PSID of the model to be
%           loaded (def: all)
% verbose : a logical if verbose should be printed (def: false)
%
% OUTPUT
% elInclMat : an binary matrix gxm model indicating whether element g in the
%               ref model is present in model m
% elRef : a cell array of g elements in the ref model
% PSID  : a cell array of m PSID that indexes the columns of geneInclMat
%
% Francesco Gatto - 2016-04-08
%

% Check args
if nargin < 1
    throw(MException('','A type of matrix must be provided among "genes","rxns","mets"'))
end
if nargin < 2 || isempty(status)
    status = 'T';
end
if nargin < 3 || isempty(PSID)
    PSID = {};
end
if nargin < 4 || isempty(verbose)
    verbose = false;
end

%Load ref model
load('../data/cHMR3765_20140415.mat')
geneRef = cModel.genes;
rxnRef  = cModel.rxns;
metRef  = cModel.mets;

%Load reconstructed models
if isempty(PSID)
    switch status
        case 'T'
            dir = '../reconstructGEMs/Tumor/Models/MAT - Copy (gatto@chalmers.se)/';
            fid = fopen(strcat(dir,'listGoodReconstructedModelMATnames.txt'),'r');
        case 'N'
            dir = '../normalGEMs/Models/';
            fid = fopen(strcat(dir,'listGoodReconstructedModelMATnames_N.txt'),'r');
    end    
    t   = textscan(fid,'%s\n');
    fclose(fid);
    modelFileNames = t{:};
    PSID=cellfun(@(a) a(1:28),modelFileNames,'uni',false);
end
nModels = numel(PSID);

%Build reference matrix
switch type
    case 'genes'
        elInclMat = zeros(numel(geneRef),nModels);
        for p =1:numel(PSID)
            models  = importModels(PSID(p),status,verbose);
            elInclMat(:,p) = ismember(geneRef,models.genes);
        end
        elRef = geneRef;
    case 'rxns'
        elInclMat = zeros(numel(rxnRef),nModels);
        for p =1:numel(PSID)
            models  = importModels(PSID(p),status,verbose);
            elInclMat(:,p) = ismember(rxnRef,models.rxns);
        end
        elRef = rxnRef;
    case 'mets'
        elInclMat = zeros(numel(metRef),nModels);
        for p =1:numel(PSID)
            models  = importModels(PSID(p),verbose);
            elInclMat(:,p) = ismember(metRef,models.mets);
        end   
        elRef = metRef;     
end