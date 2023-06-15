
# bbbPackage

### Table of Contents

1. [Installation](#installation)
2. [Functionalities](#functionalities)
3. [Examples](#examples)
4. [Files](#files)

## Installation <a name="installation"></a>
After downloading the bbbPackage-1.0-py3.whl file from the dist folder of the repository, the package can be<br/>
installed via the command pip install filepath/bbbPackage-1.0-py3.whl. Where filepath corresponds to the path <br/>
to the folder the .whl file is located in.<br/>
Additional packages used in this package will be installed as well and consist of: <br/> 
joblib, numpy, pandas, openpyxl, rdkit, mordred <br/>

## Functionalities <a name="Functionalities"></a>
This package consists of the functionalities necessary to create the Blood-Brain barrier (bbb) MI-DSE<br/> 
Fingerprints described in "PAPERNAME" for a given SMILES or set of SMILES and predict their bbb-permeability.<br/>

Function ProdFP() produces the aforementioned fingerprints and accepts either a single SMILES string, a .csv<br/> 
file containing SMILES strings and their corresponding IDs, or a .sdf file containing molecules and their<br/>
corresponding IDs. In both the latter cases activities can be supplied as well.<br/>
ProdFP() has the following keyword parameters:<br/> 
*smiles*:<br/> 
a single input SMILES string, if a set containing multiple molecules should be fingerprinted use filepath.<br/> 
By default an empty string.<br/>
*filepath*:<br/>
filepath specifying location of the .csv or .sdf file used as input. By default an empty string.<br/> 
*id_name*:<br/>
only needed in case a .sdf file is used as input. Specifies the property of the molecules holding their ID.<br/>
By default an empty string.<br/>
*act_name*:<br/>
only needed in case a .sdf file is used as input. Specifies the property of the molecules holding their activity.<br/>
By default an empty string.<br/> 
*smiles_pos*:<br/>
only needed in case a .csv file is used as input. Specifies the position of the SMILES string in the file's rows.<br/>
By default -1.<br/> 
*id_pos*:<br/>
only needed in case a .csv file is used as input. Specifies the position of the IDs in the file's rows.<br/>
By default -1.<br/> 
*act_pos*:<br/>
only needed in case a .csv file is used as input. Specifies the position of the activity in the file's rows.<br/>
By default -1.<br/> 
*skip_first*:<br/>
only needed in case a .csv file is used as input. Whether to skip first line, in case column names are present.<br/>
By default False.<br/>
ProdFP() returns a pandas Dataframe holding the descriptors the fingerprints are based on, the fingerprints<br/>
as well as the corresponding IDs and activities.<br/>

The function BbbPred() predicts bbb-permeability based on the fingerprints contained in the output Dataframe<br/>
of ProdFP(). Its parameters are<br/>
*inp_df*:<br/>
corresponds to the Dataframe returned by ProdFP().<br/>
*act*:<br/>
if set to True, BbbPred() will produce a .csv, and a .xlsx file detailing the performance of the model.<br/>
Only in case an input file containing activities was supplied can this option be used. By default False.<br/>
*ret*:<br/>
if set to True, BbbPred() returns the IDs, predictions, and prediction probabilities as lists.<br/>
Otherwise each molecule's ID, the corresponding prediction and prediction probability are printed.<br/>
By default False.<br/>
## Examples <a name="Examples"></a>
```python
import os
import bbbPackage
<!-- specifying filepath -->
filepath = os.path.join("dir1","input.csv")



```
## Files <a name="Files"></a>





















































