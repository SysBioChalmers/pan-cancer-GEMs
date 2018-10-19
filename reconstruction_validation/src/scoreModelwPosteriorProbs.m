function [rxnScores,geneScores] = scoreModelwPosteriorProbs(model,genes,posteriorProbs,noGeneScore,multipleGeneScoring)
%Remove all genes from the structures that are not in model or that aren't
%measured in the tissue
if ~isempty(genes) %This should not be necessary, but the summation is a 0x1 matrix and the other is []
    I=~ismember(genes,model.genes);
else
    I=[];
end
genes(I)=[];
posteriorProbs(I) = [];
scores = posteriorProbs;

%Get the scores for the genes, only use HPA if available
geneScores=inf(numel(model.genes),1)*-1;
[I J]=ismember(model.genes,genes);
geneScores(I)=scores(J(I));

%Remove the genes that have no data from the model
I=ismember(model.genes,genes);
model.genes(~I)=[];
model.rxnGeneMat(:,~I)=[];

%Map the genes to the HPA/array genes
[hpaExist hpaMap]=ismember(model.genes,genes);

%Set the default scores for reactions without genes
rxnScores=ones(numel(model.rxns),1)*noGeneScore;

%Loop through the reactions and calculate the scores
for i=1:numel(model.rxns)
   %Check if it has genes
   I=find(model.rxnGeneMat(i,:));
   if any(I)
       %If any of the genes exist in hpaData, then don't use arrayData
       if any(hpaExist(I))
           %At least one gene was found in HPA
           if strcmpi(multipleGeneScoring,'best')
               rxnScores(i)=max(scores(hpaMap(I(hpaExist(I)))));
           else
               rxnScores(i)=prod(scores(hpaMap(I(hpaExist(I)))));
           end
       end
   end
end
end