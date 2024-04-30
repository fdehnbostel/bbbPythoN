
# bbbPythoN

**License:** [CC BY-NC-ND 4.0](https://creativecommons.org/licenses/by-nc-nd/4.0/legalcode)

### Table of Contents

1. [Overview](#overview)
2. [Installation](#installation)
2. [Functionalities](#functionalities)
3. [Examples](#examples)
4. [Folders and Files](#Folders and Files)


## Overview <a name="overview"></a>
bbbPythoN is a Python package for the prediction of Blood-Brain barrier (BBB) permeation of small molecules.<br/>
Based on fingerprints encoding values of a pre-selected set of descriptors chosen for their discriminatory <br/>
power in regard to BBB permeability, random forest classifiers predict whether compounds are likely to <br/>
cross the BBB. Two different models, trained on an imbalanced and a balanced dataset, respectively, are <br/>
included  in this package. <br/>
As it was shown in "Non-Animal Models for Blood-Brain Barrier Permeability Evaluation of <br/>
Drug-Like Compounds" our models are most suited for the prediction of passive permeation.<br/>
Given a set of compounds, either as SMILES via command line input, as .csv file containing SMILES, or as <br/>
.sdf file containing MolBlocks, fingerprints can be produced and the molecules' activity predicted via<br/>
the random forests.<br/>

## Installation <a name="installation"></a>
After downloading the bbbPackage-1.0-py3.whl file from the dist folder of the repository, the package can be<br/>
installed via the command pip install filepath/bbbPackage-1.0-py3.whl. Where filepath corresponds to the path <br/>
to the folder the .whl file is located in.<br/>
Additional packages used in this package will be installed as well and consist of: <br/> 
joblib, numpy, pandas, openpyxl, rdkit, mordred <br/>

## Functionalities <a name="Functionalities"></a>
This package consists of the functionalities necessary to create the Blood-Brain barrier (bbb) MI-DSE<br/> 
Fingerprints described in "Non-Animal Models for Blood-Brain Barrier Permeability Evaluation of Drug-<br/>
Like Compounds" for a given SMILES, set of SMILES or set of molecules supplied via .sdf and predict<br/>
their bbb-permeability.<br/>

Function ProdFP() produces the aforementioned fingerprints and accepts either a single SMILES string, a .csv<br/> 
file containing SMILES strings and their corresponding IDs, or a .sdf file containing molecules and their<br/>
corresponding IDs. In both the latter cases activities can be supplied as well.<br/>
ProdFP() returns a pandas Dataframe holding the descriptors the fingerprints are based on, the fingerprints<br/>
as well as the corresponding IDs and activities.<br/>
ProdFP() has the following keyword parameters:<br/> 
\
*smiles*:<br/> 
A single input SMILES string, if a set containing multiple molecules should be fingerprinted use filepath.<br/> 
By default an empty string.<br/>
\
*filepath*:<br/>
Filepath specifying location of the .csv or .sdf file used as input. By default an empty string.<br/> 
\
*id_name*:<br/>
Either specifies the molecule property holding the IDs in the .sdf file.<br/>
Or the name of the column in the .csv file containing the IDs.<br/>
By default an empty string.<br/>
\
*act_name*:<br/>
Either specifies the molecule property holding the activities in the .sdf file.<br/>
Or the name of the column in the .csv file containing the activities.<br/>
By default an empty string.<br/> 
\
*smiles_name*:<br/>
Only needed in case a .csv file is used as input. Specifies the name of the <br/>
column containing the SMILES strings in the .csv file.<br/>
\
*use_bal*:<br/>
Specifies whether to use the boundaries for binning of descriptor values based <br/>
on the balanced, or imabalanced dataset.<br/>
By default False.<br/>
\
The function BbbPred() predicts bbb-permeability based on the fingerprints contained in the output Dataframe<br/>
of ProdFP(). Its parameters are:<br/>
\
*inp_df*:<br/>
Corresponds to the Dataframe returned by ProdFP().<br/>
\
*act*:<br/>
If set to True, BbbPred() will produce a .csv, and a .xlsx file detailing the performance of the model.<br/>
Only in case an input file containing activities was supplied can this option be used. By default False.<br/>
\
*ret*:<br/>
If set to True, BbbPred() returns the IDs, predictions, and prediction probabilities as lists.<br/>
Otherwise each molecule's ID, the corresponding prediction and prediction probability are printed.<br/>
By default False.<br/>

*use_bal*:<br/>
Specifies whether to use the model trained on the balanced, or imbalanced dataset.<br/>
By default False.<br/>

## Examples <a name="Examples"></a>
```python
import os
from bbbPythoN import bbbPythoN
# Defining filepath of input .csv. 
filepath_csv = os.path.join("dir1","input.csv")
# Defining filepath of input .sdf.
filepath_sdf = os.path.join("dir2","input.sdf")
# As mordred uses parallelization in its descriptor calculation the "if __name__ == '__main__':" statement
# is necessary to not spawn subprocesses from outside the main script.
if __name__ == '__main__':
# Calculating fingerprints and specifying column names of SMILES, ID, and activity in .csv.
# In case no activities are available omit act_name. 
# Exchange "id_name", "act_name", and "smiles_name" with the respective column names 
# of your dataset.
	csv_fp_df = bbbPythoN.ProdFP(filepath=filepath_csv,id_name="id_name",act_name="act_name",smiles_name="act_name",use_bal=True)
# Calculating fingerprints and molecule properties for compounds contained in .sdf file.
# In case no activities are available omit act_name.
# Exchange "id_name", "act_name", and "smiles_name" with the respective column names 
# of your dataset.
	sdf_fp_df = bbbPythoN.ProdFP(filepath=filepath_sdf,id_name="id_name",act_name="act_name",smiles_name="act_name",use_bal=True)
# It is advised to save the pandas DataFrames using .to_pickle() instead of .to_csv(), as the latter 
# will cast the python list objects representing the Fingerprints to python strings. 
	csv_fp_df.to_pickle(os.path.join("dir1","output.pkl"))
	sdf_fp_df.to_pickle(os.path.join("dir2","output.pkl"))
# Predicting bbb-permeability of .csv compounds. 
# If act=True, a .csv file listing True Positives, True Negatives, False Positives, and 
# False Negatives and a .xlsx containing performance metrics are produced.
# Setting ret=True, returns the IDs, predictions (result), and prediction probabilities (probas).
	ids, result, probas = bbbPythoN.BbbPred(csv_fp_df,act=True,ret=True,use_bal=True)
# If ret=False the molecules's IDs, predictions, and prediction probabilities are printed.
	bbbPythoN.BbbPred(csv_fp_df,ret=False,use_bal=True)
```
## Folders and Files <a name="Folders and Files"></a>

*bin_bounds*:<br/>
Contains .csv files specifying the boundaries of bins used to encode descriptor values<br/>
into bit strings for every descriptor. *bin_bounds_bal.csv*, and *bin_bounds_imba.csv* hold <br/>
boundaries based on the balanced, and the imbalanced dataset, respectively. <br/>
\
*data*:<br/>
Contains the different datasets used in "Non-Animal Models for Blood-Brain Barrier Permeability<br/>
Evaluation of Drug-Like Compounds". 
*B3DB_BalTrainTestData.csv*: the balanced dataset used to train and validate models.<br/>
*FeatSel_Data.csv*: the dataset used for MI-DSE feature selection.<br/>
*Inhouse_ImbaTrainTestData.csv*: the imbalanced dataset used to train and validate models.<br/>
*Li_Tong_Wang_ExtVal.csv*: the additonal external validation dataset.<br/>
\
*desc_names*:<br/>
Contains *mi-dse_descs.csv*, and *na_descs.csv* specifying the names of descriptors chosen by<br/>
the MI-DSE feature selection, and the names of descriptors chosen based on their discrepancy<br/>
of relative frequency of absence between classes, respectively.<br/>
\
*mah_dist_info*:<br/>
Contains *train_act_mah-dist_info_2D+3D_ma_bal*, and *train_act_mah-dist_info_2D+3D_ma_imba*,<br/>
as well as *train_mah-dist_info_2D+3D_ma_bal*, and *train_mah-dist_info_2D+3D_ma_imba*.<br/>
The former two contain information necessary for the compuation of the Mahalanobis distance <br/>
between the center of the training set *actives* and the query molecules.<br/>
The latter two contain information necessary for the compuation of the Mahalanobis distance <br/>
between the center of the training set and the query molecules.<br/>  
\
*model*:<br/>
Contains the two Random Forest classifiers trained on the balanced and imbalanced datasets.<br/>
*bbbRf_bal.sav*: model trained and validated on the balanced dataset.<br/>
*bbbRf_imba.sav*: model trained and validated on the imbalanced dataset.<br/>


















































