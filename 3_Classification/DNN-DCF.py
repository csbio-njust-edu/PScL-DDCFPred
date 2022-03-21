import argparse
import numpy as np
import pandas as pd
import sys
import pickle
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import StratifiedKFold
from sklearn.metrics import accuracy_score
from scipy.io import loadmat, matlab, savemat
import statistics
import normalize_data
sys.path.insert(0, "lib")
import utils.tools as utils

from keras.layers import Dense, merge,Input,Dropout,Flatten
from keras.models import Model
from xgboost.sklearn import XGBClassifier

from gcforest.gcforest import GCForest
from gcforest.utils.config_utils import load_json
from numpy.random import seed
seed(13)
import tensorflow
tensorflow.random.set_seed(13)

def to_categorical(y, nb_classes=None):
    '''Convert class vector (integers from 0 to nb_classes)
    to binary class matrix, for use with categorical_crossentropy.
    '''
    y = np.array(y, dtype='int')
    if not nb_classes:
        nb_classes = np.max(y)+1
    Y = np.zeros((len(y), nb_classes))
    for i in range(len(y)):
        Y[i, y[i]] = 1.
    return Y

def get_con_model(n):
    input_1 = Input(shape=(n,), name='Protein')
    protein_input1 = Dense(300, activation='relu', kernel_initializer='glorot_normal', name='High_dim_feature_1')(input_1)
    protein_input1=Dropout(0.5)(protein_input1)
    protein_input1 = Dense(200, activation='relu', kernel_initializer='glorot_normal', name='High_dim_feature_2')(protein_input1)
    protein_input1=Dropout(0.5)(protein_input1)
    output = Dense(100, activation='relu', kernel_initializer='glorot_normal', name='High_dim_feature')(protein_input1)
    outputs = Dense(6, activation='softmax', name='output')(output)
    model = Model(inputs=input_1, outputs=outputs)
    model.compile(loss='categorical_crossentropy', optimizer='Adam', metrics=['accuracy'])
    return model

def parse_args():
    parser = argparse.ArgumentParser()
    # parser.add_argument("--model", dest="model", type=str, default='E://PhD//01-RESEARCH-Bioinformatics//DCF//gcForest-master//my_gcForest.json', help="gcfoest Net Model File")
    parser.add_argument("--model", dest="model", type=str, default=None, help="gcfoest Net Model File")
    args = parser.parse_args()
    return args


def get_toy_config():
    config = {}
    ca_config = {}
    ca_config["random_state"] = 13
    ca_config["max_layers"] = 10
    ca_config["early_stopping_rounds"] = 5
    ca_config["n_classes"] = 6
    ca_config["estimators"] = []
    ca_config["estimators"].append(
            {"n_folds": 10,"booster": "gbtree", "type": "XGBClassifier", "n_estimators": 500, "max_depth": 10,
              "objective": "multi:softprob","eval_metric": "mlogloss","colsample_bytree": 0.4,"subsample": 0.8, "nthread": -1, "learning_rate": 0.1,"reg_alpha": 0.003,"gamma":0.02,
              "alpha":5,"use_label_encoder":False} )
    ca_config["estimators"].append({"n_folds": 10, "type": "RandomForestClassifier", "n_estimators": 500, "max_depth": 10, "n_jobs": -1})
    ca_config["estimators"].append({"n_folds": 10, "type": "ExtraTreesClassifier", "n_estimators": 500, "max_depth": 10, "n_jobs": -1})
    config["cascade"] = ca_config
    return config

def load_mat(filename):
    """
    This function should be called instead of direct scipy.io.loadmat
    as it cures the problem of not properly recovering python dictionaries
    from mat files. It calls the function check keys to cure all entries
    which are still mat-objects
    """

    def _check_vars(d):
        """
        Checks if entries in dictionary are mat-objects. If yes
        todict is called to change them to nested dictionaries
        """
        for key in d:
            if isinstance(d[key], matlab.mio5_params.mat_struct):
                d[key] = _todict(d[key])
            elif isinstance(d[key], np.ndarray):
                d[key] = _toarray(d[key])
        return d

    def _todict(matobj):
        """
        A recursive function which constructs from matobjects nested dictionaries
        """
        d = {}
        for strg in matobj._fieldnames:
            elem = matobj.__dict__[strg]
            if isinstance(elem, matlab.mio5_params.mat_struct):
                d[strg] = _todict(elem)
            elif isinstance(elem, np.ndarray):
                d[strg] = _toarray(elem)
            else:
                d[strg] = elem
        return d

    def _toarray(ndarray):
        """
        A recursive function which constructs ndarray from cellarrays
        (which are loaded as numpy ndarrays), recursing into the elements
        if they contain matobjects.
        """
        if ndarray.dtype != 'float64':
            elem_list = []
            for sub_elem in ndarray:
                if isinstance(sub_elem, matlab.mio5_params.mat_struct):
                    elem_list.append(_todict(sub_elem))
                elif isinstance(sub_elem, np.ndarray):
                    elem_list.append(_toarray(sub_elem))
                else:
                    elem_list.append(sub_elem)
            return np.array(elem_list)
        else:
            return ndarray

    data = loadmat(filename, struct_as_record=False, squeeze_me=True)
    return _check_vars(data)

if __name__ == "__main__":
    args = parse_args()
    if args.model is None:
        config = get_toy_config()
    else:
        config = load_json(args.model)

    gc = GCForest(config)
    dcf = GCForest(config)
    # If the model you use cost too much memory for you.
    # You can use these methods to force gcforest not keeping model in memory
    # gc.set_keep_model_in_mem(False), default is TRUE.


## Load Optimal Feats
    Data_SLFs = load_mat('sda_gda_DNAHAR.mat')
    SLFs = Data_SLFs['DNAHAR_red']
    
    Data_LBP = load_mat('sda_gda_LBP.mat')
    lbp = Data_LBP['LBP_red']
    
    Data_CLBP = load_mat('sda_gda_CLBP.mat')
    clbp = Data_CLBP['CLBP_red']
    
    
    Data_RICLBP = load_mat('sda_gda_RICLBP.mat')
    riclbp = Data_RICLBP['RICLBP_red']
    
    Data_LET = load_mat('sda_gda_LET.mat')
    let = Data_LET['LET_red']

allacc = []

class RIC_SL_CL_LE_LB_SDAGDA():
    pass
item = [RIC_SL_CL_LE_LB_SDAGDA() for j in range(10)]

i=0

for fold in SLFs:
    if (fold != 'KFold'):
        trainlabel = SLFs[fold]['trainlabel']
        testlabel = SLFs[fold]['testlabel']
        trainlabel=trainlabel-1
        testlabel=testlabel-1
        
        SLFstr = SLFs[fold]['Hartr_SDA_GDA']
        SLFste = SLFs[fold]['Harte_SDA_GDA']
        
    
        LBPtr = lbp[fold]['LBPtr_SDA_GDA']
        LBPte = lbp[fold]['LBPte_sDA_GDA']
        
        
        CLBPtr = clbp[fold]['CLBPtr_SDA_GDA']
        CLBPte = clbp[fold]['CLBPte_sDA_GDA']
        
		
        RICLBPtr = riclbp[fold]['RICLBPtr_SDA_GDA']
        RICLBPte = riclbp[fold]['RICLBPte_SDA_GDA']
        
        
        LETtr = let[fold]['LETtr_SDA_GDA']
        LETte = let[fold]['LETte_SDA_GDA']
        

        
        trainX = np.c_[SLFstr,LBPtr,CLBPtr,RICLBPtr,LETtr]
        testX = np.c_[SLFste,LBPte,CLBPte, RICLBPte,LETte]
 


        M_trainY=to_categorical(trainlabel)
        M_testY=to_categorical(testlabel)
        [m,input_dim]=np.shape(trainX)
        out_dim = 6
        cv_clf =get_con_model(input_dim)
        hist=cv_clf.fit(trainX, 
                        M_trainY,
                        epochs=100)
        y_score_train = cv_clf.predict(trainX)
        y_score_test=cv_clf.predict(testX)
        
        ##
        train_set = np.c_[trainX,y_score_train]
        test_set = np.c_[testX,y_score_test]
        
        
        X_train_enc = gc.fit_transform(train_set, trainlabel)
        
        scores = gc.predict_proba(test_set)
		y_pred = utils.categorical_probas_to_classes(scores)
        acc = accuracy_score(testlabel, y_pred)
        print("Test Accuracy of GcForest = {:.2f} %".format(acc * 100))
        allacc.append(acc)
        
        
        item[i].RI_SL_CL_LE_LB_SDAGDAte_Y = testlabel
        item[i].RI_SL_CL_LE_LB_SDAGDAteY_pred = y_pred
        item[i].RI_SL_CL_LE_LB_SDAGDAteScores = scores
        i+=1


savemat('Final_DNNDCF.mat', {'RIC_SL_CL_LE_LB_SDAGDA' : item})        
print("Mean Test Accuracy DCF:%.3f%%" %(statistics.mean(allacc)*100))