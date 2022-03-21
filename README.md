# PScL-DDCFPred
The code is the implementation of our method described in the paper “ Matee Ullah, Fazal Hadi, Jiangning Song, Dong-Jun Yu, PScL-DDCFPred: an ensemble deep learning-based approach for characterizing multiclass subcellular localization of human proteins from bioimage data”.
## (I) 1_Feature_Extraction
### (1)	data
There are two datasets:
#### (i)	Train dataset
The benchmark training dataset contains a total of 3,727 immunohistochemistry (IHC) images of six different protein subcellular locations selected from the human protein atlas (HPA) database.
#### (ii)	Independent dataset
The independent dataset contains a total of 166 IHC images of six different proteins selected from HPA. <br />
Please download the datasets from "https://drive.google.com/drive/folders/13MLIzhgcSzSCnmpsnrlyZujLPxnbZYkB?usp=sharing" and copy it to "data" folder.
### (2)	lib
lib folder contains all the code used in this study.<br />
### (3)	Biomage_Feature_Extraction.m
Biomage_Feature_Extraction.m is the matlab file for extracting <br />
(1) Subcellular location features (Slfs) which includes <br />
	(i)		DNA distribution features <br />
	(ii)	Haralick texture features <br />
(2)	Local binary pattern <br />
(3)	Completed local binary patterns <br />
(4)	Rotation invariant co-occurrence of adjacent LBP <br />
(5)	Locally encoded transform feature histogram and <br /> 
## (II)	2_Feature_Selection
2_Feature_Selection folder includes the following files
### (1)	featureSelectionCode
featureSelectionCode folder contains all the feature selection related codes and libraries.
### (2) SDA_GDA.m
SDA_GDA.m is the matlab file for the implementatoin of SDA-GDA feature selection technique. It selects the optimal features from each extracted feature set.
## (III)	3_Classification
3_Classification folder includes the following files
### (1)	lib
lib folder contains all the python codes and libraries for deep cascade forest. Complete details are available at
* Repository: https://github.com/LAMDA-NJU/Deep-Forest
* Documentation: https://deep-forest.readthedocs.io/
* Package on PyPI: https://pypi.org/project/deep-forest/
### (2)	utils
utils fodler contains the important libraries reltated to DNN-DCF classifer.
### (3)	DNN-DCF.py
DNN-DCF.py is the Python file for the implementation of DNN-DCF classifier.
## (IX)	Contact
If you are interested in our work or if you have any suggestions and questions about our research work, please contact us. E-mail: khan_bcs2010@hotmail.com.
