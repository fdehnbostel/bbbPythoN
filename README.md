
# bbbPythoN

**License:** [CC BY-NC-ND 4.0](https://creativecommons.org/licenses/by-nc-nd/4.0/legalcode)

### Table of Contents

1. [Overview](#overview)
2. [Installation](#installation)
2. [Functionalities](#functionalities)
3. [Examples](#examples)
4. [Files](#files)


## Overview <a name="overview"></a>
bbbPythoN is a Python package for the prediction of Blood-Brain barrier (BBB) permeation of small molecules.<br/>
Based on fingerprints encoding values of a pre-selected set of descriptors chosen for their discriminatory <br/>
power in regard to BBB permeability, a random forest classifier predicts whether compounds are likely to <br/>
cross the BBB. As it was shown in "Non-Animal Models for Blood-Brain Barrier Permeability Evaluation of <br/>
Drug-Like Compounds" our model is most suited for the prediction of passive permeation.<br/>
Given a set of compounds, either as SMILES via command line input, a .csv file containing SMILES, or as <br/>
.sdf file containing MolBlocks, fingerprints can be produced and the molecules' activity predicted via<br/>
the random forest.<br/>

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
a single input SMILES string, if a set containing multiple molecules should be fingerprinted use filepath.<br/> 
By default an empty string.<br/>
\
*filepath*:<br/>
filepath specifying location of the .csv or .sdf file used as input. By default an empty string.<br/> 
\
*id_name*:<br/>
only needed in case a .sdf file is used as input. Specifies the property of the molecules holding their ID.<br/>
By default an empty string.<br/>
\
*act_name*:<br/>
only needed in case a .sdf file is used as input. Specifies the property of the molecules holding their activity.<br/>
By default an empty string.<br/> 
\
*smiles_pos*:<br/>
only needed in case a .csv file is used as input. Specifies the position of the SMILES string in the file's rows.<br/>
By default -1.<br/> 
\
*id_pos*:<br/>
only needed in case a .csv file is used as input. Specifies the position of the IDs in the file's rows.<br/>
By default -1.<br/> 
\
*act_pos*:<br/>
only needed in case a .csv file is used as input. Specifies the position of the activity in the file's rows.<br/>
By default -1.<br/>
\
*skip_first*:<br/>
only needed in case a .csv file is used as input. Whether to skip first line, in case column names are present.<br/>
By default False.<br/>


\
The function BbbPred() predicts bbb-permeability based on the fingerprints contained in the output Dataframe<br/>
of ProdFP(). Its parameters are<br/>
\
*inp_df*:<br/>
corresponds to the Dataframe returned by ProdFP().<br/>
\
*act*:<br/>
if set to True, BbbPred() will produce a .csv, and a .xlsx file detailing the performance of the model.<br/>
Only in case an input file containing activities was supplied can this option be used. By default False.<br/>
\
*ret*:<br/>
if set to True, BbbPred() returns the IDs, predictions, and prediction probabilities as lists.<br/>
Otherwise each molecule's ID, the corresponding prediction and prediction probability are printed.<br/>
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
# Calculating fingerprints and specifying position of SMILES, ID, and activity in .csv rows.
# In case no activities are available omit act_pos. 
# In case a first row containing column names is present set skip_first=True.
	csv_fp_df = bbbPythoN.ProdFP(filepath=filepath_csv,smiles_pos=1,id_pos=0,act_pos= 2,skip_first=False)
# Calculating fingerprints and molecule properties for compounds contained in .sdf file.
# In case no activities are available omit act_name.
	sdf_fp_df = bbbPythoN.ProdFP(filepath=filepath_sdf,id_name="ID",act_name="Act")
# It is advised to save the pandas DataFrames using .to_pickle() instead of .to_csv(), as the latter 
# will cast the Fingerprints that are saved as python list objects to python strings. 
	csv_fp_df.to_pickle(os.path.join("dir1","output.pkl"))
	sdf_fp_df.to_pickle(os.path.join("dir2","output.pkl"))
# Predicting bbb-permeability of .csv compounds. 
# If act=True, a .csv file listing True Positives, True Negatives, False Positives, and 
# False Negatives and a .xlsx containing performance metrics are produced.
# Setting ret=True, returns the IDs, predictions (result), and prediction probabilities (probas).
	ids, result, probas = bbbPythoN.BbbPred(csv_fp_df,act=True,ret=True)
# If ret=False the molecules's IDs, predictions, and prediction probabilities are printed.
	bbbPythoN.BbbPred(csv_fp_df,ret=False)
```
## Files <a name="Files"></a>

bin_bounds.csv: contains the boundaries of bins used to encode descriptor values into bit strings<br/>
for every descriptor.<br/>
\
mi-dse_descs.csv: contains the names of the descriptors chosen by MI-DSE feature selection.<br/>
\
na_descs.csv: contains the names of descriptors chosen based on their discrepancy in relative<br/> 
frequency of absence between classes.<br/>
\
train_act_mah-dist_info_2D+3D_ma: contains information necessary for the compuation of the<br/> 
Mahalanobis distance between the center of the training set actives and the query molecules.<br/> 
\
train_mah-dist_info_2D+3D_ma: contains information necessary for the compuation of the<br/> 
Mahalanobis distance between the center of the training set and the query molecules.<br/> 


















































