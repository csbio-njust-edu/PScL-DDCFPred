clear;
clc;

load all_labels;
load all_features_SLFs_LBPs;
load all_features_CLBP;
load all_features_LET
load all_features_RICLBP

% labels_TR = label_all;
SLFs = all_features_SLFs_LBPs(:,2:841);
LBP = all_features_SLFs_LBPs(:,842:end);
CLBP = all_features_CLBP;
RICLBP = all_features_RICLBP;
LET = all_features_LET;
SLFs(isnan(SLFs)|isinf(SLFs))=0;
LBP(isnan(LBP)|isinf(LBP))=0;
CLBP(isnan(CLBP)|isinf(CLBP))=0;
RICLBP(isnan(RICLBP)|isinf(RICLBP))=0;
LET(isnan(LET)|isinf(LET))=0;
CVO = cvpartition(label_all', 'k', 10);

SLFsoptions.KernelType = 'rbf';
SLFsoptions.KernelPars = 4.3;

LBPoptions.KernelType = 'rbf';
LBPoptions.KernelPars = 1.3;

CLBPoptions.KernelType = 'rbf';
CLBPoptions.KernelPars = 2.3;

RICLBPoptions.KernelType = 'rbf';
RICLBPoptions.KernelPars = 1.7;

LEToptions.KernelType = 'rbf';
LEToptions.KernelPars = 3.9;

for fold=1:10
    trIdxUn=[]; teIdxUn=[]; DNAtr=[]; DNAte=[]; SLFstr=[]; SLFste=[];
    LBPtr=[]; LBPte=[]; CLBPtr=[]; CLBPte=[]; trainlabel=[]; testlabel=[]; 
    TRsim = []; TEsim = []; DecisionValuetr = []; DecisionValue = [];
    
    trIdxUn = CVO.training(fold);
    teIdxUn = CVO.test(fold);

    kfold = strcat('kfold',num2str(fold));
    
    SLFstr = SLFs(trIdxUn,:);
    SLFste = SLFs(teIdxUn,:);
    LBPtr = LBP(trIdxUn,:);
    LBPte = LBP(teIdxUn,:);
    CLBPtr = CLBP(trIdxUn,:);
    CLBPte = CLBP(teIdxUn,:);
	LETtr = LET(trIdxUn,:);
    LETte = LET(teIdxUn,:);
    RICLBPtr = RICLBP(trIdxUn,:);
    RICLBPte = RICLBP(teIdxUn,:);
    trainlabel= label_all(trIdxUn,:);
    testlabel =label_all(teIdxUn,:);
    
    
%% SLFs Normalization
    massimo=max(SLFstr);
    minimo=min(SLFstr);
    for i=1:size(SLFstr,2)
        SLFstr(1:size(SLFstr,1),i)=double(SLFstr(1:size(SLFstr,1),i)-minimo(i))/(massimo(i)-minimo(i));
    end
    for i=1:size(SLFste,2)
        SLFste(1:size(SLFste,1),i)=double(SLFste(1:size(SLFste,1),i)-minimo(i))/(massimo(i)-minimo(i));
    end
    
%% LBP Normalization
    massimo=max(LBPtr);
    minimo=min(LBPtr);
    for i=1:size(LBPtr,2)
        LBPtr(1:size(LBPtr,1),i)=double(LBPtr(1:size(LBPtr,1),i)-minimo(i))/(massimo(i)-minimo(i));
    end
    for i=1:size(LBPte,2)
        LBPte(1:size(LBPte,1),i)=double(LBPte(1:size(LBPte,1),i)-minimo(i))/(massimo(i)-minimo(i));
    end
    
%% CLBP Normalization
    massimo=max(CLBPtr);
    minimo=min(CLBPtr);
    for i=1:size(CLBPtr,2)
        CLBPtr(1:size(CLBPtr,1),i)=double(CLBPtr(1:size(CLBPtr,1),i)-minimo(i))/(massimo(i)-minimo(i));
    end
    for i=1:size(CLBPte,2)
        CLBPte(1:size(CLBPte,1),i)=double(CLBPte(1:size(CLBPte,1),i)-minimo(i))/(massimo(i)-minimo(i));
    end
	
	% LET Normalization
    massimo=max(LETtr);
    minimo=min(LETtr);
    for i=1:size(LETtr,2)
        LETtr(1:size(LETtr,1),i)=double(LETtr(1:size(LETtr,1),i)-minimo(i))/(massimo(i)-minimo(i));
    end
    for i=1:size(LETte,2)
        LETte(1:size(LETte,1),i)=double(LETte(1:size(LETte,1),i)-minimo(i))/(massimo(i)-minimo(i));
    end
    
% RICLBP Normalization
    massimo=max(RICLBPtr);
    minimo=min(RICLBPtr);
    for i=1:size(RICLBPtr,2)
        RICLBPtr(1:size(RICLBPtr,1),i)=double(RICLBPtr(1:size(RICLBPtr,1),i)-minimo(i))/(massimo(i)-minimo(i));
    end
    for i=1:size(RICLBPte,2)
        RICLBPte(1:size(RICLBPte,1),i)=double(RICLBPte(1:size(RICLBPte,1),i)-minimo(i))/(massimo(i)-minimo(i));
    end
 
%% SLFs Feature selection
    [idx_sdaSLFs SLFstr SLFste] = SDA_FeatSelect(double(SLFstr), double(SLFste), trainlabel);
    
    DNASLFs_red.(kfold).idx_sda = idx_sdaSLFs;
    DNASLFs_red.(kfold).SLFstr_SDA = SLFstr;
    DNASLFs_red.(kfold).SLFste_sDA = SLFste;
    DNASLFs_red.(kfold).trainlabel = trainlabel;
    DNASLFs_red.(kfold).testlabel = testlabel;
    
        SLFstr = double(SLFstr');
        SLFste = double(SLFste');
        trainlabel = trainlabel';
        SLFstrgda  =  gda(SLFstr,SLFstr,trainlabel,5,SLFsoptions);     % Project the training data matrix into a low-dimensional space
        SLFstegda  =  gda(SLFste,SLFstr,trainlabel,5,SLFsoptions);       % Project the test data matrix into a low-dimensional space
        SLFstr = SLFstrgda';
        SLFste = SLFstegda';
        trainlabel = trainlabel';
    DNASLFs_red.(kfold).SLFstr_SDA_GDA = SLFstr;
    DNASLFs_red.(kfold).SLFste_SDA_GDA = SLFste;
 
%% LBP Feature selection 
     [idx_sdaLBP LBPtr LBPte] = SDA_FeatSelect(LBPtr, LBPte, trainlabel);
    
    LBP_red.(kfold).idx_sda = idx_sdaLBP;
    LBP_red.(kfold).LBPtr_SDA = LBPtr;
    LBP_red.(kfold).LBPte_sDA = LBPte;
    LBP_red.(kfold).trainlabel = trainlabel;
    LBP_red.(kfold).testlabel = testlabel;
    
        LBPtr = double(LBPtr');
        LBPte = double(LBPte');
        trainlabel = trainlabel';
        LBPtrgda  =  gda(LBPtr,LBPtr,trainlabel,5,LBPoptions);     % Project the training data matrix into a low-dimensional space
        LBPtegda  =  gda(LBPte,LBPtr,trainlabel,5,LBPoptions);       % Project the test data matrix into a low-dimensional space
        LBPtr = LBPtrgda';
        LBPte = LBPtegda';
        trainlabel = trainlabel';
    
    LBP_red.(kfold).LBPtr_SDA_GDA = LBPtr;
    LBP_red.(kfold).LBPte_sDA_GDA = LBPte;
    
%% CLBP Feature selection
     [idx_sdaCLBP CLBPtr CLBPte] = SDA_FeatSelect(double(CLBPtr), double(CLBPte), trainlabel);
     
    CLBP_red.(kfold).idx_sda = idx_sdaCLBP;
    CLBP_red.(kfold).CLBPtr_SDA = CLBPtr;
    CLBP_red.(kfold).CLBPte_sDA = CLBPte;
    CLBP_red.(kfold).trainlabel = trainlabel;
    CLBP_red.(kfold).testlabel = testlabel;
     
        CLBPtr = double(CLBPtr');
        CLBPte = double(CLBPte');
        trainlabel = trainlabel';
        CLBPtrgda  =  gda(CLBPtr,CLBPtr,trainlabel,5,CLBPoptions);     % Project the training data matrix into a low-dimensional space
        CLBPtegda  =  gda(CLBPte,CLBPtr,trainlabel,5,CLBPoptions);       % Project the test data matrix into a low-dimensional space
        CLBPtr = CLBPtrgda';
        CLBPte = CLBPtegda';
        trainlabel = trainlabel';
        
    CLBP_red.(kfold).CLBPtr_SDA_GDA = CLBPtr;
    CLBP_red.(kfold).CLBPte_sDA_GDA = CLBPte;

%% LET Feature selection
[idx_sdaLET LETtr LETte] = SDA_FeatSelect(double(LETtr), double(LETte), trainlabel);

    LET_red.(kfold).idx_sda = idx_sdaLET;
    LET_red.(kfold).LETtr_SDA = LETtr;
    LET_red.(kfold).LETte_sDA = LETte;
    LET_red.(kfold).trainlabel = trainlabel;
    LET_red.(kfold).testlabel = testlabel;


        LETtr = LETtr';
        LETte = LETte';
        trainlabel = trainlabel';
        LETtrGda  =  gda(LETtr,LETtr,trainlabel,5,LEToptions);     % Project the training data matrix into a low-dimensional space
        LETteGda  =  gda(LETte,LETtr,trainlabel,5,LEToptions);       % Project the test data matrix into a low-dimensional space
        trainlabel = trainlabel';
        
    LET_red.(kfold).LETtr_SDA_GDA = LETtrGda';
    LET_red.(kfold).LETte_SDA_GDA = LETteGda';
    
%% RICLBP Feature selection
[idx_sdaRIC RICLBPtr RICLBPte] = SDA_FeatSelect(double(RICLBPtr), double(RICLBPte), trainlabel);

    RICLBP_red.(kfold).idx_sda = idx_sdaRIC;
    RICLBP_red.(kfold).RICLBPtr_SDA = RICLBPtr;
    RICLBP_red.(kfold).RICLBPte_sDA = RICLBPte;
    RICLBP_red.(kfold).trainlabel = trainlabel;
    RICLBP_red.(kfold).testlabel = testlabel;

	
        RICLBPtr = RICLBPtr';
        RICLBPte = RICLBPte';
        trainlabel = trainlabel';
        RICLBtrGda  =  gda(RICLBPtr,RICLBPtr,trainlabel,5,RICLBPoptions);     % Project the training data matrix into a low-dimensional space
        RICLBteGda  =  gda(RICLBPte,RICLBPtr,trainlabel,5,RICLBPoptions);       % Project the test data matrix into a low-dimensional space
        trainlabel = trainlabel';
        
    RICLBP_red.(kfold).RICLBPtr_SDA_GDA = RICLBtrGda';
    RICLBP_red.(kfold).RICLBPte_SDA_GDA = RICLBteGda';	
end
save sda_gda_SLFs SLFs_red;
save sda_gda_LBP LBP_red;
save sda_gda_CLBP CLBP_red;
save sda_gda_LET LET_red;
save sda_gda_RICLBP RICLBP_red;