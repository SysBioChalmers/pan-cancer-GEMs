%% Generate model from post probs
function TINITmodel = getINITmodelFromPostProbs(patientIndex,probabilityThreshold,relativeGap)
% getINITmodelFromPostProbs
%
% This function generates a tINIT model from PostProbs in a file. refModel, 
% scores and tasks must be specified within this script. 

%%% Initialize path
ravenPath  = genpath('/path/to/RAVEN/');
mosek7path = genpath('/path/to/mosek/7');
addpath(ravenPath);
addpath(mosek7path);

% Load data
load ../cHMR3765_20140415
load ../allFiles_AllMutCID_BSfactors_minDE190.mat
load ../cHMR3765_Tasks     %Task structure
if exist('hotStart.mat','file')==2
    load ../hotStart.mat                   %MILP Hot Start
    hotStart              = res.sol.int.xx;
else
    hotStart = [];
end

%%% Files
refModel     = cModel;
description  = strcat(PSID{patientIndex},'T',num2str(probabilityThreshold*10000),'RG',num2str(relativeGap*100));
presentMets  = [];
if exist('../cHMR3765_essentialRxnsForTasks.mat','file')==2
    load ../cHMR3765_essentialRxnsForTasks.mat
else
    load ../cHMR3765_Tasks
    essentialRxnsForTasks = getEssentialRxnsForTasks(refModel, [], taskList, true);
end
taskStructure = taskList;

%%% Algorithm settings
allowExcretion      = true;  %Enforce presence of reactions that may be regarded as dead-ends just because a complete biomass eqn is missing
printReport         = true;
params.MSK_DPAR_OPTIMIZER_MAX_TIME = 96*60*60;
params.MSK_DPAR_MIO_TOL_REL_GAP    = relativeGap;
paramsFT.maxTime      = 96*60;
paramsFT.printReport  = true;
paramsFT.relGap       = relativeGap;

%%% Score settings
posteriorProbs      = pP(patientIndex,:) - probabilityThreshold;
genes               = ENSG;
noGeneScore         = 0;
multipleGeneScoring = 'best';
useScoresForTasks   = true;
allowNewMets        = false;

%%% Run INIT and tINIT
fprintf('Running INIT...\n')
tStartINIT = tic;
INITmodel = getINITModelInServer(refModel, description, presentMets,...
    printReport, params, allowExcretion, essentialRxnsForTasks,hotStart,...
    posteriorProbs,genes,noGeneScore,multipleGeneScoring);
fprintf('\tINIT successfully terminated.\n')
tElapsedINIT=toc(tStartINIT);
fprintf('Running tINIT...\n')
tStartTINIT = tic;
[TINITmodel addedRxnsForTasks] = getTINITModelInServer(refModel, INITmodel, description, ...
    useScoresForTasks, printReport, taskStructure, paramsFT, allowNewMets, ...
    posteriorProbs,genes,noGeneScore,multipleGeneScoring);
fprintf('\tINIT successfully terminated.\n')
tElapsedTINIT=toc(tStartTINIT);

%%% Evaluate model
[~,deletedDeadEndRxns]     = simplifyModel(refModel,true,false,true,true,true);
nodeadendModel             = removeRxns(refModel,deletedDeadEndRxns,false,true);
[~,geneRefScores]          = scoreModelwPosteriorProbs(nodeadendModel,genes,pP(patientIndex,:),noGeneScore,multipleGeneScoring);
geneScoresInINIT           = geneRefScores(ismember(nodeadendModel.genes,TINITmodel.genes));
geneScoresNotInINIT        = geneRefScores(~ismember(nodeadendModel.genes,TINITmodel.genes));
TP  = sum(geneScoresInINIT > probabilityThreshold);
FP  = sum(geneScoresInINIT <= probabilityThreshold);
TN  = sum(geneScoresNotInINIT <= probabilityThreshold);
FN  = sum(geneScoresNotInINIT > probabilityThreshold);
MCC = (TP*TN - FP*FN)/(sqrt((TP+FP)*(TP+FN)*(TN+FP)*(TN+FN)));
ACC = (TP+TN)/(TP+TN+FP+FN);
p   = chi2pdf(MCC^2*(TP+TN+FP+FN),1);
percAddedRxnsForTasks = numel(addedRxnsForTasks)/numel(TINITmodel.rxns)*100;

%%% Save
outputfile = strcat(description,'_TINIT.mat');
    save(outputfile,'TINITmodel')
    exportModel(TINITmodel,strcat(TINITmodel.id,'.xml'),true) %Save as .xml COBRA
fid = fopen('INITmodelSummary.txt','a');
fprintf(fid,'%s\t%s\t%s\t%d\t%f\t%f\t%f\t%f\t%e\t%f\t%f\t%f\n',outputfile,strcat(TINITmodel.id,'.xml'),...
    PSID{patientIndex},patientIndex,probabilityThreshold,relativeGap,ACC,MCC,p,percAddedRxnsForTasks,tElapsedINIT,tElapsedTINIT);
fclose(fid);
fprintf('Saved.\n')

end

function [initModel metProduction deletedDeadEndRxns deletedRxnsInINIT]=getINITModelInServer(refModel, description, ...
    presentMets, printReport, params, allowExcretion, essentialRxnsForTasks,hotStart,...
    posteriorProbs,genes,noGeneScore,multipleGeneScoring)
% getINITModelInServer
%   Generates a model using the INIT algorithm, based on proteomics and/or
%   transcriptomics and/or metabolomics and/or metabolic tasks.
%
%   refModel            a model structure. The model should be in the
%                       closed form (no exchange reactions open). Import
%                       using import(filename,false). If the model is not
%                       loaded using importModel, it might be that there
%                       is no "unconstrained" field. In that case,
%                       manually add the field like:
%                       model.unconstrained=false(numel(model.mets),1);
%   description         a string describing the resulting model.
%   presentMets         cell array with metabolite names that the model
%                       should produce (opt, default [])
%   printReport         true if a report should be printed to the screen
%                       (opt, default true)
%   params              parameter structure as used by getMILPParams. This is
%                       for the INIT algorithm. For the the MILP problems
%                       solved to fit tasks, see paramsFT (opt, default [])
%   allowExcretion    true if excretion of all metabolites should be allowed.
%                     This results in fewer reactions being considered
%                     dead-ends, but all reactions in the resulting model may
%                     not be able to carry flux. If this is "false" then the
%                     equality constraints are taken from model.b. If the
%                     input model lacks exchange reactions then this should
%                     probably be "true", or a large proportion of the model
%                     would be excluded for connectivity reasons
%                     (opt, default false)
%   essentialRxnsForTasks   a cell array of rxn IDs that are essential to
%                           perform a list of predefined tasks (opt, default [])
%   hotStart                a MOSEK structure containing a hot-start for
%                           MILP optimization (opt, default [])
%   posteriorProbs          a numeric array nx1 containing the posterior
%                           probabilities of n genes to be expressed in the
%                           given condition
%   genes                   a cell array of n genes that index the
%                           posterior probabilities and should be contained
%                           in refModel.genes
%   noGeneScore             a numerical indicating the score for genes in
%                           refModel.genes that are not accounted for in
%                           posteriorProbs
%   multipleGeneScoring     a string either 'best' or 'average' that
%                           indicate how to treat posterior probabilities
%                           for multiple genes refering to the same rxn.
%
%   model                   the resulting model structure
%   metProduction           array that indicates which of the
%                           metabolites in presentMets that could be
%                           produced. Note that this is before the
%                           gap-filling process to enable defined tasks. To
%                           see which metabolites that can be produced in
%                           the final model, use canProduce. 
%                           -2: metabolite name not found in model
%                           -1: metabolite found, but it could not be produced
%                           1: metabolite could be produced
%   deletedDeadEndRxns      cell array of reactions deleted because they
%                           could not carry flux (INIT requires a
%                           functional input model)
%   deletedRxnsInINIT       cell array of the reactions which were deleted by 
%                           the INIT algorithm

%   This is the main function for automatic reconstruction of models based
%   on the INIT algorithm (PLoS Comput Biol. 2012;8(5):e1002518). Not all
%   settings are possible using this function, and you may want to call the
%   functions scoreModel, runINIT and fitTasksInServer individually instead.
%
%   NOTE: Exchange metabolites should normally not be removed from the model
%   when using this approach, since checkTasks/fitTasksInServer rely on putting specific
%   constraints for each task. The INIT algorithm will remove exchange metabolites
%   if any are present. Use importModel(file,false) to import a model with
%   exchange metabolites remaining.
%
%   Usage: [initModel metProduction deletedDeadEndRxns deletedRxnsInINIT]=getINITModelInServer(refModel, description, ...
%    presentMets, printReport, params, allowExcretion, essentialRxnsForTasks,hotStart,...
%    posteriorProbs,genes,noGeneScore,multipleGeneScoring)
%   
%   Rasmus Agren, 2013-05-16
%   Francesco Gatto, 2014-06-12

%Check that the model is in the closed form
if ~isfield(refModel,'unconstrained')
   throw(MException('','Exchange metabolites should normally not be removed from the model when using getINITModel. Use importModel(file,false) to import a model with exchange metabolites remaining (see the documentation for details)')); 
end

if printReport==true
        fprintf(['***Generating model for: ' description '\n']);
    if ~isempty(presentMets)
        fprintf('-Using metabolomics data\n');
    end
    fprintf('\n');
    printScores(refModel,'Reference model statistics',genes,posteriorProbs,noGeneScore,multipleGeneScoring);
end


%Remove dead-end reactions to speed up the optimization and to
%differentiate between reactions removed by INIT and those that are
%dead-end
[crap, deletedDeadEndRxns] = simplifyModel(refModel,true,false,true,true,true);
cModel=removeRxns(refModel,deletedDeadEndRxns,false,true);

%Store the connected model like this to keep track of stuff
if printReport==true
    printScores(cModel,'Pruned model statistics',genes,posteriorProbs,noGeneScore,multipleGeneScoring);
end

%Score the connected model
[rxnScores,geneScores] = scoreModelwPosteriorProbs(cModel,genes,posteriorProbs,noGeneScore,multipleGeneScoring);

%Run the INIT algorithm. The exchange reactions that are used in the final
%reactions will be open, which doesn't fit with the last step. Therefore
%delete reactions from the original model instead of taking the output.
%The default implementation does not constrain reversible reactions to only
%carry flux in one direction.
%Runs without the constraints on reversibility and with all output allowed.
%This is to reduce the complexity of the problem.

[crap deletedRxnsInINIT metProduction]=runINITinServer(simplifyModel(cModel),rxnScores,presentMets,essentialRxnsForTasks,0,allowExcretion,false,params,hotStart);
initModel = removeRxns(cModel,deletedRxnsInINIT,true,true);
if printReport==true
    printScores(initModel,'INIT model statistics',genes,posteriorProbs,noGeneScore,multipleGeneScoring);
    printScores(removeRxns(cModel,setdiff(cModel.rxns,deletedRxnsInINIT),true,true),'Reactions deleted by INIT',genes,posteriorProbs,noGeneScore,multipleGeneScoring);
end
%The full model has exchange reactions in it. fitTasksInServer calls on fillGaps,
%which automatically removes exchange metabolites (because it assumes that
%the reactions are constrained when appropriate). In this case the
%uptakes/outputs are retrieved from the task sheet instead. To prevent
%exchange reactions being used to fill gaps, they are delete from the
%reference model here.
initModel.id= description;
end

%This is for printing a summary of a model
function [rxnS geneS]=printScores(model,name,genes,posteriorProbs,noGeneScore,multipleGeneScoring)
    [a b]= scoreModelwPosteriorProbs(model,genes,posteriorProbs,noGeneScore,multipleGeneScoring);
    rxnS=mean(a);
    geneS=mean(b(~isinf(b)));
    fprintf([name ':\n']);
    fprintf(['\t' num2str(numel(model.rxns)) ' reactions, ' num2str(numel(model.genes)) ' genes\n']);
    fprintf(['\tMean reaction score: ' num2str(rxnS) '\n']);
    fprintf(['\tMean gene score: ' num2str(geneS) '\n']);
    fprintf(['\tReactions with positive scores: ' num2str(100*sum(a>0)/numel(a)) '%%\n\n']);
end

%This is to save the mosek problem and give a hot start
function [outModel deletedRxns metProduction fValue]=runINITinServer(model,rxnScores,presentMets,essentialRxns,prodWeight,allowExcretion,noRevLoops,params,hotStart)
% runINIT
%	Generates a model using the INIT algorithm, based on proteomics and/or
%   transcriptomics and/or metabolomics and/or metabolic tasks
%
%   model           a reference model structure
%   rxnScores       a vector of scores for the reactions in the model.
%                   Positive scores are reactions to keep and negative
%                   scores are reactions to exclude (opt, default all 0.0)
%   presentMets     cell array with unique metabolite names that the model
%                   should produce (opt, default [])
%   essentialRxns   cell array of reactions that are essential and that
%                   have to be in the resulting model. This is normally
%                   used when fitting a model to task (see fitTasks) (opt,
%                   default [])
%   prodWeight      a score that determines the value of having
%                   net-production of metabolites. This is a way of having
%                   a more functional network as it provides a reason for
%                   including bad reactions for connectivity reasons. This
%                   score is for each metabolite, and the sum of these weights
%                   and the scores for the reactions is what is optimized
%                   (opt, default 0.5)
%   allowExcretion  true if excretion of all metabolites should be allowed.
%                   This results in fewer reactions being considered
%                   dead-ends, but all reactions in the resulting model may
%                   not be able to carry flux. If this is "false" then the
%                   equality constraints are taken from model.b. If the
%                   input model lacks exchange reactions then this should
%                   probably be "true", or a large proportion of the model
%                   would be excluded for connectivity reasons
%                   (opt, default false)
%   noRevLoops      true if reversible reactions should be constrained to
%                   only carry flux in one direction. This prevents
%                   reversible reactions from being wrongly assigned as
%                   connected (the forward and backward reactions can form a
%                   loop and therefore appear connected), but it makes the 
%                   problem significantly more computationally intensive to
%                   solve (two more integer constraints per reversible reaction)
%                   (opt, default false)
%   params          parameter structure as used by getMILPParams (opt,
%                   default [])
%   hotStart        an array containing a feasible integer solution to the
%                   problem
%
%   outModel        the resulting model structure
%   deletedRxns     reactions which were deleted by the algorithm
%   metProduction   array that indicates which of the
%                   metabolites in presentMets that could be
%                   produced 
%                   -2: metabolite name not found in model
%                   -1: metabolite found, but it could not be produced
%                   1: metabolite could be produced
%   fValue          objective value (sum of (the negative of) 
%                   reaction scores for the included reactions and
%                   prodWeight*number of produced metabolites)
%
%   This function is the actual implementation of the algorithm. See
%   getINITModel for a higher-level function for model reconstruction. See 
%   PLoS Comput Biol. 2012;8(5):e1002518 for details regarding the
%   implementation.
%
%   Usage: [outModel deletedRxns metProduction fValue]=runINITinServer(model,...
%           rxnScores,presentMets,essentialRxns,prodWeight,allowExcretion,...
%           noRevLoops,params,hotStart)
%
%   Rasmus Agren, 2013-08-01
%

if nargin<2
    rxnScores=zeros(numel(model.rxns),1);
end
if isempty(rxnScores)
    rxnScores=zeros(numel(model.rxns),1);
end
if nargin<3
    presentMets={};
end
if isempty(presentMets)
    presentMets={};
end
presentMets=presentMets(:);
if nargin<4
    essentialRxns={};
end
if isempty(essentialRxns)
    essentialRxns={};
end
essentialRxns=essentialRxns(:);
if nargin<5
    prodWeight=0.5;
end
if isempty(prodWeight)
    prodWeight=0.5;
end
if nargin<6
    allowExcretion=false;
end
if nargin<7
    noRevLoops=false;
end
if nargin<8
    params=[];
end
if nargin<9
    hotStart = [];
end

echo=0;
if isfield(params,'printReport')
    if params.printReport==true
        echo=3;
    end
end

if numel(presentMets)~=numel(unique(presentMets))
	dispEM('Duplicate metabolite names in presentMets');
end

%Default is that the metabolites cannot be produced
if ~isempty(presentMets)
    metProduction=ones(numel(presentMets),1)*-2;
    presentMets=upper(presentMets);
    pmIndexes=find(ismember(presentMets,upper(model.metNames)));
    metProduction(pmIndexes)=-1; %Then set that they are at least found
else
    metProduction=[];
    pmIndexes=[];
end

%The model should be in the reversible format and all relevant exchange
%reactions should be open
if isfield(model,'unconstrained')
    dispEM('Exchange metabolites are still present in the model. Use simplifyModel if this is not intended',false); 
end

%The irreversible reactions that are essential must have a flux and are therefore not
%optimized for using MILP, which reduces the problem size. However, reversible
%reactions must have a flux in one direction, so they have to stay in
%the problem. The essentiality constraint on reversible reactions is
%implemented in the same manner as for reversible reactions when
%noRevLoops==true, but with the additional constraint that C ub=-1. This
%forces one of the directions to be active.
revRxns=find(model.rev~=0);
essentialReversible=find(ismember(model.rxns(revRxns),essentialRxns));
essentialRxns=intersect(essentialRxns,model.rxns(model.rev==0));

%Convert the model to irreversible
irrevModel=convertToIrrev(model);
rxnScores=[rxnScores;rxnScores(model.rev==1)];
%These are used if noRevLoops is true
if noRevLoops==true
    forwardIndexes=find(model.rev~=0);
    backwardIndexes=(numel(model.rxns)+1:numel(irrevModel.rxns))';
else
    %Then they should only be used for essential reversible reactions
    forwardIndexes=revRxns(essentialReversible);
    backwardIndexes=essentialReversible+numel(model.rxns);
end

%Get the indexes of the essential reactions and remove them from the
%scoring vector
essentialIndex=find(ismember(irrevModel.rxns,essentialRxns));
rxnScores(essentialIndex)=[];

%Go through each of the presentMets (if they exist) and modify the S matrix
%so that each reaction which produces any of them also produces a
%corresponding fake metabolite and the opposite in the reverse direction.

%This is to deal with the fact that there is no compartment info regarding
%the presentMets. This modifies the irrevModel structure, but that is fine
%since it's the model structure that is returned.
if any(pmIndexes)
    irrevModel.metNames=upper(irrevModel.metNames);
    metsToAdd.mets=strcat({'FAKEFORPM'},num2str(pmIndexes));
    metsToAdd.metNames=metsToAdd.mets;
    metsToAdd.compartments=irrevModel.comps{1};
    
    %There is no constraints on the metabolites yet, since maybe not all of
    %them could be produced
    irrevModel=addMets(irrevModel,metsToAdd);
end

%Modify the matrix
for i=1:numel(pmIndexes)
    %Get the matching mets
    I=ismember(irrevModel.metNames,presentMets(pmIndexes(i)));
    
    %Find the reactions where any of them are used.
    [crap K L]=find(irrevModel.S(I,:));
    
    %This ugly loop is to avoid problems if a metabolite occurs several
    %times in one reaction
    KK=unique(K);
    LL=zeros(numel(KK),1);
    for j=1:numel(KK)
       LL(j)=sum(L(K==KK(j)));
    end
    irrevModel.S(numel(irrevModel.mets)-numel(pmIndexes)+i,KK)=LL;
end

%Some nice to have numbers
nMets=numel(irrevModel.mets);
nRxns=numel(irrevModel.rxns);
nEssential=numel(essentialIndex);
nNonEssential=nRxns-nEssential;
nonEssentialIndex=setdiff(1:nRxns,essentialIndex);
S=irrevModel.S;

%Add so that each non-essential reaction produces one unit of a fake metabolite
temp=sparse(1:nRxns,1:nRxns,1);
temp(essentialIndex,:)=[];
S=[S;temp];

%Add another set of reactions (will be binary) which also produce these
%fake metabolites, but with a stoichiometry of 1000
temp=sparse(1:nNonEssential,1:nNonEssential,1000);
temp=[sparse(nMets,nNonEssential);temp];
S=[S temp];

%Add reactions for net-production of (real) metabolites
if prodWeight~=0
    temp=[speye(nMets-numel(pmIndexes))*-1;sparse(nNonEssential+numel(pmIndexes),nMets-numel(pmIndexes))];
    S=[S temp];
    %To keep the number of reactions added like this
    nNetProd=nMets-numel(pmIndexes);
else
    nNetProd=0;
end

%Add constraints so that reversible reactions can only be used in one
%direction. This is done by adding the fake metabolites A, B, C for each
%reversible reaction in the following manner
% forward: A + .. => ...
% backwards: B + ... => ...
% int1: C => 1000 A
% int2: C => 1000 B
% A ub=999.9
% B ub=999.9
% C lb=-1
% int1 and int2 are binary
if any(forwardIndexes)
    nRevBounds=numel(forwardIndexes);
    
    %Add the A metabolites for the forward reactions and the B
    %metabolites for the reverse reactions
    I=speye(numel(irrevModel.rxns))*-1;
    temp=[I(forwardIndexes,:);I(backwardIndexes,:)];
   
    %Padding
    temp=[temp sparse(size(temp,1),size(S,2)-numel(irrevModel.rxns))];
    
    %Add the int1 & int2 reactions that produce A and B
    temp=[temp speye(nRevBounds*2)*1000];
    
    %And add that they also consume C
    temp=[temp;[sparse(nRevBounds,size(S,2)) speye(nRevBounds)*-1 speye(nRevBounds)*-1]];
    
    %Add the new reactions and metabolites
    S=[S sparse(size(S,1),nRevBounds*2)];
    S=[S;temp];
else
    nRevBounds=0;
end

%Add so that the essential reactions must have a small flux and that the
%binary ones (and net-production reactions) may have zero flux. The
%integer reactions for reversible reactions have [0 1]
prob.blx=[irrevModel.lb;zeros(nNonEssential+nNetProd+nRevBounds*2,1)];
prob.blx(essentialIndex)=max(0.1,prob.blx(essentialIndex));

%Add so that the binary ones and net-production reactions can have at the most flux 1.0
prob.bux=[irrevModel.ub;ones(nNonEssential+nNetProd+nRevBounds*2,1)];

%Add that the fake metabolites must be produced in a small amount and that
%the A and B metabolites for reversible reactions can be [0 999.9] and C
%metabolites [-1 0]
prob.blc=[irrevModel.b(:,1);ones(nNonEssential,1);zeros(nRevBounds*2,1);ones(nRevBounds,1)*-1];

%Add that normal metabolites can be freely excreted if allowExcretion==true,
%and that the fake ones can be excreted 1000 units at most. C metabolites
%for essential reversible reactions should have an upper bound of -1.
%If noRevLoops is false, then add this constraint for all the reactions instead.
if noRevLoops==true
    revUB=zeros(nRevBounds,1);
    revUB(essentialReversible)=-1;
else
    revUB=ones(nRevBounds,1)*-1;
end
if allowExcretion==true
    metUB=inf(nMets,1);
else
    metUB=irrevModel.b(:,min(size(irrevModel.b,2),2));
end
prob.buc=[metUB;ones(nNonEssential,1)*1000;ones(nRevBounds*2,1)*999.9;revUB];

%Add objective coefficients for the binary reactions. The negative is used
%since we're minimizing. The negative is taken for the prodWeight
%as well, in order to be consistent with the syntax that positive scores
%are good
prob.c=[zeros(nRxns,1);rxnScores;ones(nNetProd,1)*prodWeight*-1;zeros(nRevBounds*2,1)];
prob.a=S;

%We still don't know which of the presentMets that can be produced. Go
%through them, force production, and see if the problem can be solved
%%% WARNING %%%
%%% params.MSK_IPAR_OPTIMIZER='MSK_OPTIMIZER_FREE_SIMPLEX'; - TRY OUT %%%
%%% DEFAULT SOLVER %%%
for i=1:numel(pmIndexes)
    prob.blc(numel(irrevModel.mets)-numel(pmIndexes)+i)=1;
    [crap,res] = mosekopt('minimize echo(0)', prob,params);
    isFeasible=checkSolution(res);
    if ~isFeasible
        %Reset the constraint again
        prob.blc(numel(irrevModel.mets)-numel(pmIndexes)+i)=0;
    else
        %Metabolite produced
        metProduction(pmIndexes(i))=1;
    end
end

%Add that the binary reactions may only take integer values.
allInt=[(nRxns+1):(nRxns+nNonEssential) size(S,2)-nRevBounds*2+1:size(S,2)];
prob.ints.sub=allInt;

%Get params and provide a hot start if available
% mosekParam = getMILPParams(params);
save('prob.mat','prob');
if ~isempty(hotStart)
    if size(prob.a,2) ~= max(size(hotStart))
        fprintf('WARNING: hotStart must have the same dimensions as the MILP formulation. The hot start will be ignored\n')
    else
        prob.sol.int.xx = hotStart;
    end
end
[crap,res] = mosekopt(['minimize echo(' num2str(echo) ')'], prob,params);

%I don't think that this problem can be infeasible, so this is mainly a way
%of checking the licence stuff
[isFeasible isOptimal] = checkSolution(res);
if ~isFeasible
    dispEM('The problem is infeasible');
end
if ~isOptimal
    dispEM('The problem is not optimal');
end

fValue=res.sol.int.pobjval;

%Get all reactions used in the irreversible model
usedRxns=(nonEssentialIndex(res.sol.int.xx(nRxns+1:nRxns+nNonEssential)<0.1))';

%Map to reversible model IDs
usedRxns=[usedRxns(usedRxns<=numel(model.rxns));revRxns(usedRxns(usedRxns>numel(model.rxns))-numel(model.rxns))];

%Then get the ones that are not used in either direction or is essential
I=true(numel(model.rxns),1);
I(usedRxns)=false;
I(essentialIndex)=false;
deletedRxns=model.rxns(I);
outModel=removeRxns(model,I,true,true);
end

%This is to partition the computations carried out in mosek7 or 6 and to
%allow for scoring with array
function [model,addedRxnsForTasks] = getTINITModelInServer(refModel, initModel, description, ...
    useScoresForTasks, printReport, taskStructure, paramsFT, allowNewMets, ...
    posteriorProbs,genes,noGeneScore,multipleGeneScoring)

% getTINITModelInServer
%   Generates a model using the INIT algorithm, based on proteomics and/or
%   transcriptomics and/or metabolomics and/or metabolic tasks.
%
%   refModel            a model structure. The model should be in the
%                       closed form (no exchange reactions open). Import
%                       using import(filename,false). If the model is not
%                       loaded using importModel, it might be that there
%                       is no "unconstrained" field. In that case,
%                       manually add the field like:
%                       model.unconstrained=false(numel(model.mets),1);
%   tissue              tissue to score for. Should exist in either
%                       hpaData.tissues or arrayData.tissues
%   celltype            cell type to score for. Should exist in either
%                       hpaData.celltypes or arrayData.celltypes for this 
%                       tissue (opt, default is to use the best values
%                       among all the cell types for the tissue. Use [] if
%                       you want to supply more arguments)
%   hpaData             HPA data structure from parseHPA (opt if arrayData is 
%                       supplied, default [])
%   arrayData           gene expression data structure (opt if hpaData is
%                       supplied, default [])
%       genes           cell array with the unique gene names
%       tissues         cell array with the tissue names. The list may not be
%                       unique, as there can be multiple cell types per tissue
%       celltypes       cell array with the cell type names for each tissue
%       levels          GENESxTISSUES array with the expression level for
%                       each gene in each tissue/celltype. NaN should be
%                       used when no measurement was performed
%   metabolomicsData    cell array with metabolite names that the model
%                       should produce (opt, default [])
%   taskFile            a task list in Excel format. See parseTaskList for
%                       details (opt, default [])
%   useScoresForTasks   true if the calculated reaction scored should be used as
%                       weights in the fitting to tasks (opt, default true)
%   printReport         true if a report should be printed to the screen
%                       (opt, default true)
%   taskStructure       task structure as from parseTaskList. Can be used
%                       as an alternative way to define tasks when Excel
%                       sheets are not suitable. Overrides taskFile (opt,
%                       default [])
%   params              parameter structure as used by getMILPParams. This is
%                       for the INIT algorithm. For the the MILP problems
%                       solved to fit tasks, see paramsFT (opt, default [])
%   paramsFT            parameter structure as used by getMILPParams. This is
%                       for the fitTasksInServer step. For the INIT algorithm, see
%                       params (opt, default [])
%   allowExcretion  true if excretion of all metabolites should be allowed.
%                     This results in fewer reactions being considered
%                     dead-ends, but all reactions in the resulting model may
%                     not be able to carry flux. If this is "false" then the
%                     equality constraints are taken from model.b. If the
%                     input model lacks exchange reactions then this should
%                     probably be "true", or a large proportion of the model
%                     would be excluded for connectivity reasons
%                     (opt, default false)
%
%   model                   the resulting model structure
%
%   This is the main function for automatic reconstruction of models based
%   on the INIT algorithm (PLoS Comput Biol. 2012;8(5):e1002518). Not all
%   settings are possible using this function, and you may want to call the
%   functions scoreModel, runINIT and fitTasksInServer individually instead.
%
%   NOTE: Exchange metabolites should normally not be removed from the model
%   when using this approach, since checkTasks/fitTasksInServer rely on putting specific
%   constraints for each task. The INIT algorithm will remove exchange metabolites
%   if any are present. Use importModel(file,false) to import a model with
%   exchange metabolites remaining.
%
%   Usage: [model metProduction essentialRxnsForTasks addedRxnsForTasks...
%               deletedDeadEndRxns deletedRxnsInINIT taskReport]=...
%               getINITModel(refModel, tissue, celltype, hpaData, arrayData,...
%               metabolomicsData, taskFile, useScoresForTasks, printReport,...
%               taskStructure, params, paramsFT)
%
%   Rasmus Agren, 2013-05-16
%


%Check that the model is in the closed form
if ~isfield(refModel,'unconstrained')
   throw(MException('','Exchange metabolites should normally not be removed from the model when using getINITModel. Use importModel(file,false) to import a model with exchange metabolites remaining (see the documentation for details)')); 
end

%Create the task structure if not supplied
if isempty(taskStructure)
   throw(MException('','Tasks must be provided as a structure.')); 
end

%Remove dead-end reactions to speed up the optimization and to
%differentiate between reactions removed by INIT and those that are
%dead-end
[crap, deletedDeadEndRxns]=simplifyModel(refModel,true,false,true,true,true);
cModel=removeRxns(refModel,deletedDeadEndRxns,false,true);

%If gaps in the model should be filled using a task list
if ~isempty(taskStructure)    
    %Remove exchange reactions and reactions already included in the INIT
    %model
    refModelNoExc=removeRxns(refModel,union(initModel.rxns,getExchangeRxns(refModel)),true,true);
    
    %At this stage the model is fully connected and most of the genes with
    %good scores should have been included. The final gap-filling should
    %take the scores of the genes into account, so that "rather bad"
    %reactions are preferred to "very bad" reactions. However, reactions
    %with positive scores will be included even if they are not connected
    %in the current formulation. Therefore, such reactions will have to be
    %assigned a small negative score instead.
    if useScoresForTasks==true
        refRxnScores = scoreModelwPosteriorProbs(cModel,genes,posteriorProbs,noGeneScore,multipleGeneScoring);
        [outModel addedRxnMat]=fitTasksInServer(initModel,refModelNoExc,[],true,min(refRxnScores,-0.1),taskStructure,paramsFT,allowNewMets);
    else
        [outModel addedRxnMat]=fitTasksInServer(initModel,refModelNoExc,[],true,[],taskStructure,paramsFT,allowNewMets);
    end
    if printReport==true
        printScores(outModel,'Functional model statistics',genes,posteriorProbs,noGeneScore,multipleGeneScoring);
        printScores(removeRxns(outModel,intersect(outModel.rxns,initModel.rxns),true,true),'Reactions added to perform the tasks',genes,posteriorProbs,noGeneScore,multipleGeneScoring);
    end
    
    addedRxnsForTasks=refModelNoExc.rxns(any(addedRxnMat,2));
else
    outModel=initModel;
    addedRxnsForTasks={};
end

%The model can now perform all the tasks defined in the task list. The
%algorithm cannot deal with gene-complexes at the moment. It is therefore
%ok to remove bad genes from a reaction (as long as at least one gene is
%kept)
model=outModel;
[~,geneScoresRef] = scoreModelwPosteriorProbs(refModel,genes,posteriorProbs,noGeneScore,multipleGeneScoring);

[crap I]=ismember(model.genes,refModel.genes); %All should be found
%This is a little weird way to make sure that only one bad gene is included
%if there are no good ones (since all -Inf==max(-Inf))
geneScoresRef(isinf(geneScoresRef))=-1000+rand(sum(isinf(geneScoresRef)),1);

model.grRules(:)={''};
for i=1:numel(model.rxns)
    ids=find(model.rxnGeneMat(i,:));
    if numel(ids)>1
       scores=geneScoresRef(I(ids));
       %Only keep the positive ones if possible
       model.rxnGeneMat(i,ids(~(scores>0 | scores==max(scores))))=0;
    end
    %Rewrite the grRules to be only OR
    if isfield(model,'grRules')
       J=find(model.rxnGeneMat(i,:));
       for j=1:numel(J)
           model.grRules{i}=[model.grRules{i} '(' model.genes{J(j)} ')'];
           if j<numel(J)
               model.grRules{i}=[model.grRules{i} ' or '];
           end
       end
    end
end

%Find all genes that are not used and delete them
I=sum(model.rxnGeneMat)==0;
model.genes(I)=[];
model.rxnGeneMat(:,I)=[];
if isfield(model,'geneShortNames')
    model.geneShortNames(I)=[];
end
if isfield(model,'geneMiriams')
    model.geneMiriams(I)=[];
end
if isfield(model,'geneFrom')
    model.geneFrom(I)=[];
end

%At this stage the model will contain some exchange reactions but probably not all
%(and maybe zero). This can be inconvenient, so all exchange reactions from the 
%reference model are added, except for those which involve metabolites that
%are not in the model.

%First delete and included exchange reactions in order to prevent the order
%from changing
model=removeRxns(model,getExchangeRxns(model));

%Create a model with only the exchange reactions in refModel
excModel=removeRxns(refModel,setdiff(refModel.rxns,getExchangeRxns(refModel)),true,true);

%Find the metabolites there which are not exchange metabolites and which do
%not exist in the output model
I=~ismember(excModel.mets,model.mets) & excModel.unconstrained==0;

%Then find those reactions and delete them
[crap J]=find(excModel.S(I,:));
excModel=removeRxns(excModel,J,true,true);

%Merge with the output model
model=mergeModels({model;excModel});
model.id= description;
model.description=['Automatically generated model for ' description];

if printReport==true
    printScores(model,'Final model statistics',genes,posteriorProbs,noGeneScore,multipleGeneScoring); 
end
end

%This is to save after each task
function [outModel addedRxns]=fitTasksInServer(model,refModel,inputFile,printOutput,rxnScores,taskStructure,params,allowNewMets)
% fitTasksInServer
%   Fills gaps in a model by including reactions from a reference model,
%   so that the resulting model can perform all the tasks in a task list
%
%   model           model structure
%   refModel        reference model from which to include reactions 
%   inputFile       a task list in Excel format. See the function
%                   parseTaskList for details (opt if taskStructure is
%                   supplied)
%   printOutput     true if the results of the test should be displayed
%                   (opt, default true)
%   rxnScores       scores for each of the reactions in the reference
%                   model. Only negative scores are allowed. The solver will
%                   try to maximize the sum of the scores for the included
%                   reactions (opt, default is -1 for all reactions)
%   taskStructure   structure with the tasks, as from parseTaskList. If
%                   this is supplied then inputFile is ignored (opt)
%   params          parameter structure as used by getMILPParams (opt)
%
%   outModel        model structure with reactions added to perform the
%                   tasks
%   addedRxns       MxN matrix with the added reactions (M) from refModel 
%                   for each task (N). An element is true if the corresponding
%                   reaction is added in the corresponding task.
%                   Failed tasks and SHOULD FAIL tasks are ignored
%
%   This function fills gaps in a model by using a reference model, so
%   that the resulting model can perform a list of metabolic tasks. The
%   gap-filling is done in a task-by-task manner, rather than solving for
%   all tasks at once. This means that the order of the tasks could influence
%   the result.
%
%   Usage: [outModel addedRxns]=fitTasksInServer(model,refModel,inputFile,printOutput,...
%           rxnScores,taskStructure,params)
%
%   Rasmus Agren, 2013-04-22
%

if nargin<4
    printOutput=true;
end
if nargin<5
    rxnScores=ones(numel(refModel.rxns),1)*-1;
end
if isempty(rxnScores)
    rxnScores=ones(numel(refModel.rxns),1)*-1;
end
if nargin<6
    taskStructure=[];
end
if nargin<7
    params=[];
end

if strcmpi(model.id,refModel.id)
	fprintf('NOTE: The model and reference model have the same IDs. The ID for the reference model was set to "refModel" in order to keep track of the origin of reactions.\n'); 
    refModel.id='refModel';
end

if any(rxnScores>=0)
	throw(MException('','Only negative values are allowed in rxnScores'));
end

%Prepare the input models a little
model.b=zeros(numel(model.mets),2);
modelMets=upper(strcat(model.metNames,'[',model.comps(model.metComps),']'));
%This is the mets in the reference model. Used if the tasks involve
%metabolites that doesn't exist in the model
largeModelMets=upper(strcat(refModel.metNames,'[',refModel.comps(refModel.metComps),']'));

if ~isfield(model,'unconstrained')
	fprintf('WARNING: Exchange metabolites should normally not be removed from the model when using checkTasks. Inputs and outputs are defined in the task file instead. Use importModel(file,false) to import a model with exchange metabolites remaining.\n'); 
end

if isempty(taskStructure)
   taskStructure=parseTaskList(inputFile); 
end

tModel=model;
addedRxns=false(numel(refModel.rxns),numel(taskStructure));
supressWarnings=false;
nAdded=0;
for i=1:numel(taskStructure)
    if ~taskStructure(i).shouldFail
        %Set the inputs
        if ~isempty(taskStructure(i).inputs)
            [I J]=ismember(upper(taskStructure(i).inputs),modelMets);
            K=ismember(upper(taskStructure(i).inputs),'ALLMETS');
            L=~cellfun('isempty',strfind(upper(taskStructure(i).inputs),'ALLMETSIN'));
            %Check that all metabolites are either real metabolites or
            %ALLMETS/ALLMETSIN
            goodMets=I|K|L;
            if ~all(goodMets)
                %Not all of the inputs could be found in the small model. Check
                %if they exist in the large model
                [found metMatch]=ismember(upper(taskStructure(i).inputs(~goodMets)),largeModelMets);
                if ~all(found)
                    throw(MException('',['Could not find all inputs in "[' taskStructure(i).id '] ' taskStructure(i).description '" in either model']));
                else
                   %Otherwise add them to the model
                   met.metNames=refModel.metNames(metMatch);
                   met.compartments=refModel.comps(refModel.metComps(metMatch));
                   
                   %Add the metabolite both to the base model and the model
                   %used in the current task
                   model=addMets(model,met);
                   tModel=addMets(tModel,met);
                   modelMets=[modelMets;upper(taskStructure(i).inputs(~goodMets))];
                end
                
                %By now the indexes might be getting a bit confusing, but
                %this is to update the indexes of the "real" metabolites to
                %point to the newly added ones
                I(~goodMets)=true; %All the bad ones are fixed at this stage
                J(~goodMets)=numel(modelMets)-numel(metMatch)+1:numel(modelMets);
            end
            if numel(J(I))~=numel(unique(J(I)))
                throw(MException('',['The constraints on some input(s) in "[' taskStructure(i).id '] ' taskStructure(i).description '" are defined more than one time']));  
            end
            %If all metabolites should be added
            if any(K)
                %Check if ALLMETS is the first metabolite. Otherwise print a
                %warning since it will write over any other constraints that
                %are set
                if K(1)==0
                    fprintf(['WARNING: ALLMETS is used as an input in "[' taskStructure(i).id '] ' taskStructure(i).description '" but it it not the first metabolite in the list. Constraints defined for the metabolites before it will be over-written\n']);
                end
                %Use the first match of ALLMETS. There should only be one, but
                %still..
                tModel.b(:,1)=taskStructure(i).UBin(find(K,1))*-1;
            end
            %If metabolites in a specific compartment should be used
            if any(L)
                L=find(L);
                for j=1:numel(L)
                    %The compartment defined
                    compartment=upper(taskStructure(i).inputs{L(j)}(11:end-1));
                    %Check if it exists in the model
                    C=find(ismember(upper(model.comps),compartment));
                    if any(C)
                        %Match to metabolites
                        tModel.b(model.metComps==C,1)=taskStructure(i).UBin(L(j))*-1;
                    else
                        throw(MException('',['The compartment defined for ALLMETSIN in "[' taskStructure(i).id '] ' taskStructure(i).description '" does not exist']));  
                    end
                end
            end
            %Then add the normal constraints
            if any(J(I))
                tModel.b(J(I),1)=taskStructure(i).UBin(I)*-1;
                tModel.b(J(I),2)=taskStructure(i).LBin(I)*-1;
            end
        end
        %Set the outputs
        if ~isempty(taskStructure(i).outputs)
            [I J]=ismember(upper(taskStructure(i).outputs),modelMets);
            K=ismember(upper(taskStructure(i).outputs),'ALLMETS');
            L=~cellfun('isempty',strfind(upper(taskStructure(i).outputs),'ALLMETSIN'));
            %Check that all metabolites are either real metabolites or
            %ALLMETS/ALLMETSIN
            goodMets=I|K|L;
            if ~all(goodMets)
                %Not all of the outputs could be found in the small model. Check
                %if they exist in the large model
                [found metMatch]=ismember(upper(taskStructure(i).outputs(~goodMets)),largeModelMets);
                if ~all(found)
                    throw(MException('',['Could not find all outputs in "[' taskStructure(i).id '] ' taskStructure(i).description '" in either model']));
                else
                   %Otherwise add them to the model
                   met.metNames=refModel.metNames(metMatch);
                   met.compartments=refModel.comps(refModel.metComps(metMatch));
                   
                   %Add the metabolite both to the base model and the model
                   %used in the current task
                   model=addMets(model,met);
                   tModel=addMets(tModel,met);
                   modelMets=[modelMets;upper(taskStructure(i).outputs(~goodMets))];
                end
                
                %By now the indexes might be getting a bit confusing, but
                %this is to update the indexes of the "real" metabolites to
                %point to the newly added ones
                I(~goodMets)=true; %All the bad ones are fixed at this stage
                J(~goodMets)=numel(modelMets)-numel(metMatch)+1:numel(modelMets);
            end
            if numel(J(I))~=numel(unique(J(I)))
                throw(MException('',['The constraints on some output(s) in "[' taskStructure(i).id '] ' taskStructure(i).description '" are defined more than one time']));  
            end
            %If all metabolites should be added
            if any(K)
                %Check if ALLMETS is the first metabolite. Otherwise print a
                %warning since it will write over any other constraints that
                %are set
                if K(1)==0
                    fprintf(['WARNING: ALLMETS is used as an output in "[' taskStructure(i).id '] ' taskStructure(i).description '" but it it not the first metabolite in the list. Constraints defined for the metabolites before it will be over-written\n']);
                end
                %Use the first match of ALLMETS. There should only be one, but
                %still..
                tModel.b(:,2)=taskStructure(i).UBout(find(K,1));
            end
            %If metabolites in a specific compartment should be used
            if any(L)
                L=find(L);
                for j=1:numel(L)
                    %The compartment defined
                    compartment=upper(taskStructure(i).outputs{L(j)}(11:end-1));
                    %Check if it exists in the model
                    C=find(ismember(upper(model.comps),compartment));
                    if any(C)
                        %Match to metabolites
                        tModel.b(model.metComps==C,2)=taskStructure(i).UBout(L(j));
                    else
                        throw(MException('',['The compartment defined for ALLMETSIN in "[' taskStructure(i).id '] ' taskStructure(i).description '" does not exist']));  
                    end
                end
            end
            %Then add the normal constraints
            if any(J(I))
                tModel.b(J(I),1)=taskStructure(i).LBout(I);
                tModel.b(J(I),2)=taskStructure(i).UBout(I);
            end
        end 
        
        %Add new rxns
        if ~isempty(taskStructure(i).equations)
            rxn.equations=taskStructure(i).equations;
            rxn.lb=taskStructure(i).LBequ;
            rxn.ub=taskStructure(i).UBequ;
            rxn.rxns=strcat({'TEMPORARY_'},num2str((1:numel(taskStructure(i).equations))'));
            tModel=addRxns(tModel,rxn,3,[],allowNewMets);
        end
        %Add changed bounds
        if ~isempty(taskStructure(i).changed)
           tModel=setParam(tModel,'lb',taskStructure(i).changed,taskStructure(i).LBrxn);
           tModel=setParam(tModel,'ub',taskStructure(i).changed,taskStructure(i).UBrxn);
        end

        %Solve and print. Display a warning if the problem is not solveable
        sol=solveLP(tModel);
        if isempty(sol.x)
            %Only do gap-filling if it cannot be solved
            failed=false;
            try
                [crap crap newRxns newModel exitFlag]=fillGaps(tModel,refModel,false,true,supressWarnings,rxnScores,params);
                if exitFlag==-2
                    fprintf(['WARNING: "[' taskStructure(i).id '] ' taskStructure(i).description '" was aborted before reaching optimality. Consider increasing params.maxTime\n\n']);
                end
            catch
                fprintf(['WARNING: "[' taskStructure(i).id '] ' taskStructure(i).description '" could not be performed for any set of reactions\n\n']);
                failed=true;
            end
            if failed==false
                if ~isempty(newRxns)
                    nAdded=nAdded+numel(newRxns);
                    
                    %Add the reactions to the base model. It is not correct to use newModel
                    %directly, as it may contain reactions/constraints that are specific to
                    %this task
                    model=mergeModels({model,removeRxns(newModel,setdiff(newModel.rxns,newRxns),true,true)},true);
                    
                    %Keep track of the added reactions
                    addedRxns(ismember(refModel.rxns,newRxns),i)=true;
                end
                if printOutput==true
                    fprintf(['[' taskStructure(i).id '] ' taskStructure(i).description ': Added ' num2str(numel(newRxns)) ' reaction(s), ' num2str(nAdded) ' reactions added in total\n']);
                end
            end
        else
            if printOutput==true
                fprintf(['[' taskStructure(i).id '] ' taskStructure(i).description ': Added 0 reaction(s), ' num2str(nAdded) ' reactions added in total\n']);
            end
        end
        supressWarnings=true;
        
        %Print the output if chosen
        if taskStructure(i).printFluxes && printOutput
            if ~isempty(sol.x)
                sol=solveLP(tModel,1);
                printFluxes(tModel,sol.x,false,10^-5,[],'%rxnID (%eqn):%flux\n');
                fprintf('\n');
            else
                %If the problem wasn't solveable then the gap-filled model
                %should be used
                if failed==false
                    sol=solveLP(newModel,1);
                    printFluxes(newModel,sol.x,false,10^-5,[],'%rxnID (%eqn):%flux\n');
                    fprintf('\n');
                end
            end
        end
        
        tModel=model;
        %Since new mets are added by merging the new reactions and not only
        %from the task sheet
        modelMets=upper(strcat(model.metNames,'[',model.comps(model.metComps),']'));
    else
        fprintf(['WARNING: "[' taskStructure(i).id '] ' taskStructure(i).description '" is set as SHOULD FAIL. Such tasks cannot be modelled using this approach and the task is therefore ignored\n\n']);
    end
end
outModel=model;
end