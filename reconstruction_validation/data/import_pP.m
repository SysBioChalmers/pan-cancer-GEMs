filename = 'pP_validation.csv';
pP       = csvread(filename,1,1);
nGenes   = size(pP,2);
nPts     = size(pP,1);
dimnames = textscan(fopen(filename),'%s','Delimiter',',');
ENSG     = dimnames{1}(2:(nGenes+1));
PSID     = dimnames{1}(strncmp('TCGA',dimnames{1},4));
%PSID     = cellfun(@(x) x(1:12),PSID_long,'UniformOutput',false)
if ~length(ENSG)==size(pP,2)
   error('Error. N (%d) genes is different than n cols in pP matrix.',length(ENSG))
end


if ~length(PSID)==size(pP,1)
   error('Error. N (%d) samples is different than n rows in pP matrix.',length(ENSG))
end
save('pP.mat','pP','ENSG','PSID')