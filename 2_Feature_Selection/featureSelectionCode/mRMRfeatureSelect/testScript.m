% test script
clear; clc;
currentFolder = pwd;
addpath(genpath(currentFolder))
dataPath    = '/feature-selection-mRMR/matlab/exampleHDR/';
dataSet     = 'mfeat.mat';
dataFile    = [dataPath, dataSet];
load(dataSet );
dataX        = rawData;
dataC       = kron((1:10)',ones(200,1));
preDataMtd  = 'binarize';
% ===== Test for candidateFeature =====
% load('/Users/kylechen/Dropbox/Github/feature-selection-mRMR/dataset/mfeat/binarizedData.mat');
% dataX = binarizedData;
% dataC = kron((1:10)',ones(200,1));
nMRMR = 10;
wrapper = 'back';
classifier = 'NB';
errThres = 1e-1;
tic;
candiFea = candidateFeature(dataX, dataC, nMRMR, classifier, errThres);
% cmptFea = compactWrapper(dataX, dataC, candiFea, classifier, wrapper);
[ fCmptFea, fErrRcd ] = compactWrapper(dataX, dataC, candiFea, classifier, 'for');
[ bCmptFea, bErrRcd ] = compactWrapper(dataX, dataC, candiFea, classifier, 'back');
toc;

hPlot = plot(fErrRcd, '-ob');
hold on;
plot(bErrRcd, '->r');
% title('', 'fontsize', 16);
xlabel('feature number', 'FontSize', 16);
ylabel('eror rate', 'FontSize', 16);
legendStr   = {'forward'; 'backward'};
legend(legendStr);
filename = strcat('rho_quantile', filenameStr{seqM,model} , '.eps');
print(h,filename,'-depsc2','-r300');




% kFold = 10;
% cvErrEst(dataX(:,candiFea), dataC, classifier, kFold)


% ===== Test for candidateFeature =====
% load fisheriris
% dataX = meas;
% dataC = [ones(50,1); 2*ones(50,1); 3*ones(50,1)];
% nMRMR = 100;
% wrapper = 'back';
% classifier = 'NB';
% candiFea = candidateFeature(dataX, dataC, nMRMR, wrapper, classifier);

% ===== Test for mRMR =====
% file = '/Users/kylechen/Documents/MATLAB/hcPengMRMR/mRMR/mRMR_0.9_compiled/test_lung_s3.csv';
% dataset = load(file);
% 
% dataC = dataset(:,1);
% dataX = dataset(:,2:end);
% 
% nSelect = 10;
% tic;
% myFea   = mRMR(dataX, dataC, nSelect);
% toc;
% tic;
% feaPeng = mrmr_mid_d(dataX, dataC, nSelect);
% toc;
% Lia = ismember(myFea,feaPeng)