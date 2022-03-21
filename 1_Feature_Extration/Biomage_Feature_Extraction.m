function Biomage_Feature_Extraction()

close all, clear all, clc;
currentFolder = pwd;
addpath(genpath(currentFolder))
task = {'parse','separate','featcalc','classification','prediction','results'};
%addPath;

antibody_ids = {'2384','4125','29804','2321','36986','3901',...
                '2645','7912','38677','30827','35593','11929','35866',...
                '1524','2018','41690','4480','26533','50453','63636','53891','54422',...
                '28136','43912','11993','17964','23099','19025','22012','19684', ...
               };

classlabels{1,1}={'Cytoplasm'};
classlabels{1,2}={'Cytoplasm'};
classlabels{1,3}={'Cytoplasm'};
classlabels{1,4}={'Cytoplasm'};
classlabels{1,5}={'Cytoplasm'};
classlabels{1,6}={'ER'};
classlabels{1,7}={'ER'};
classlabels{1,8}={'ER'};
classlabels{1,9}={'ER'};
classlabels{1,10}={'ER'};
classlabels{1,11}={'Golgi_apparatus'};
classlabels{1,12}={'Golgi_apparatus'};
classlabels{1,13}={'Golgi_apparatus'};
classlabels{1,14}={'Golgi_apparatus'};
classlabels{1,15}={'Golgi_apparatus'};
classlabels{1,16}={'Mitochondria'};
classlabels{1,17}={'Mitochondria'};
classlabels{1,18}={'Mitochondria'};
classlabels{1,19}={'Mitochondria'};
classlabels{1,20}={'Mitochondria'};
classlabels{1,21}={'Nucleus'};
classlabels{1,22}={'Nucleus'};
classlabels{1,23}={'Nucleus'};
classlabels{1,24}={'Nucleus'};
classlabels{1,25}={'Nucleus'};
classlabels{1,26}={'Vesicles'};
classlabels{1,27}={'Vesicles'};
classlabels{1,28}={'Vesicles'};
classlabels{1,29}={'Vesicles'};
classlabels{1,30}={'Vesicles'};
 

umet = {'lin','nmf'};  
dbs = {'db1'};
featType = {'SLFs_LBPs','LET','RICLBP','CLBP'};
classifyMethod = {'BR','CC'};

%% 1. parameter setting and Image statistics calculation
 % This is the code for parameters setting and
 % calculates a few basic image statistics 
 % for the purpose of identifying badly stained images.
for i1=1:length(antibody_ids)
    options.general.ANTIBODY_IDS = antibody_ids(i1);
    params = setOptions( options);
    clslbls = params.general.CLASSLABELS;
    data = prepareData( 'parse',params);
    handleData( data, params);
end
%%
    label_all=[];
     save all_labels  label_all 
     number=0;
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%% 2. Image Seperation
% This block of code is used to seperate the protein and DNA images.
id =2;
if (~strcmp(task{id},'results'))  
    for m = 1:length(featType) 
        options.featCalc.FEATTYPE = featType{m}; 
        for k=1:1    
            options.separate.UMETHOD = umet{k};        
            if strcmp(options.featCalc.FEATTYPE,'SLFs_LBPs') && strcmp(options.separate.UMETHOD,'nmf')
               continue;
            end
            for j=1:length(dbs)
               options.featCalc.FEATSET = dbs{j};
               switch task{id}
                    case 'separate'
                        for i=1:9
                            options.general.ANTIBODY_IDS = antibody_ids(i);
                            options.general.CLASSLABELS = classlabels(i);
                            params = setOptions( options);
                            number=i;
                            save numbers number
                            data = prepareData( 'separate',params);
                            handleData( data, params);
                        end                  
                  end
               end           
        end
    end
  
end     
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 3. Feature extraction     
% This block of code is used to extract features from the images.    
id =3;
options.general.CLASSLABELS = clslbls;
if (~strcmp(task{id},'results'))  
    for m = 1:length(featType)
        featName = featType{m};
        save_features(featName);
        options.featCalc.FEATTYPE = featType{m}; 
        for k=1:1    
            options.separate.UMETHOD = umet{k};        
            if strcmp(options.featCalc.FEATTYPE,'SLFs_LBPs') && strcmp(options.separate.UMETHOD,'nmf')
               continue;
            end
            for j=1:length(dbs)
               options.featCalc.FEATSET = dbs{j};
               switch task{id}
                    
                    case 'featcalc'
                        for i=1:length(antibody_ids)
                            options.general.ANTIBODY_IDS = antibody_ids(i);
                            params = setOptions( options); 
                            number=i;
                            save numbers number
                            data = prepareData( 'featcalc',params);
                            handleData( data, params);
                        end                                                             
                    otherwise
               end
                if j>0
                    break
                end
            end

        end
    end
    
end
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
% function addPath( dataroot)
% %  addPath( DATAROOT)  initializes folder paths 
% %    Input DATAROOT is by default './data'
% 
if ~exist( 'dataroot','var')
    dataroot = 'data';
end
% 
directoryFun( '.', dataroot);
directoryFun( dataroot, 'dataset');
directoryFun( dataroot, '2_separatedImages');
directoryFun( dataroot, '3_features');
% 
function directoryFun( prnt, chld)
d = [prnt '/' chld];
if ~exist( d, 'dir')
    mkdir( prnt, chld);
end
    


