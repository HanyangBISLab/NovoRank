# NovoRank
NovoRank is a machine learning/deep learning-based algorithm that post-processes the De novo sequencing results.

- NovoRank is implemented and tested with

Python==3.9 \
requirements.txt

and

DeepLC ( https://github.com/compomics/DeepLC ) \
MS-Cluster ( http://proteomics.ucsd.edu/software-tools/ms-clusterarchives ) \
CometX.exe ( this is in-house software modified to calculate XCorr only. The implementation is based on the Comet software )

All data used in the experiment can be downloaded from the 'sample' folder in the NovoRank GitHub repository and https://drive.google.com/drive/folders/13jir3-QLcAVGtUe84Tfp5utctnHRI-5Q?usp=share_link.

## Quick start for potential reviewers

A user can download pre-trained model and test sample data at https://drive.google.com/drive/folders/1ICmLBRBhJGdImi4aPwlQVS6pfKT-Z8sI?usp=share_link and run below command line for quick test:

- Create and activate an Anaconda virtual environment.
```c
conda create -n [NAME] python==3.9
conda activate [NAME]
```

<br/>

- Install all packages in the requirements.txt file with the following command.
```c
pip install -r requirements.txt
```

<br/>

- Run run_novorank.py.
```c
python run_novorank.py .\test\config_for_reviewer.txt
```

## How to use NovoRank

To use NovoRank for your datasets, you HAVE TO train your own model fitting to your datasets.

### Step 1. Preparation datasets
As an initial step, a user MUST make their datasets to fit NovoRank input standard.

#### De novo search result
Once you perform de novo search using any tools such as PEAKS, pNovo3 and DeepNovo, you MUST convert the result to below form:

Source File|Scan number|Peptide|Score
---|---|---|---|
Hela_1.mgf|10|HKPSVK|85|

Note that each column is separated by comma (comma-separated value format (CSV)).

#### Database search result
NovoRank generates positive and negative labels based on database search result from the same MS/MS spectra used in the de novo search. Therefore, it only needs for training. If a user uses pre-trained model, this file is not needed for the further step. After conducting database search, only reliable PSMs are prepared as below format:

Source File|Scan number|GT
---|---|---|
Hela_1.mgf|3|KPVGAAK| 

Note that each column is separated by comma (comma-separated value format (CSV)).

#### Note for post-translational modification notation
NovoRank assumes that all Cysteines (C) have a fixed modification Carbamidomethylation.
As a variable modification, it only allows an oxidation on Methionine as lower letter "m".
For example, if AM+15.99EENGR, a user must convert the sequence to AmEENGR.

### Step 2. Initial clustering using MS-Cluster

MS-Cluster software and user’s manual are available at http://proteomics.ucsd.edu/software-tools/ms-clusterarchives/. Create a list of the full paths to the input files and call it list.txt. 

< Clustering to MS-Cluster using the following command line. >
```c
MSCluster.exe --list list.txt --output-name CLUSTERS --mixture-prob 0.01 --fragment-tolerance 0.02 --assign-charges
``` 

### Step 3. Generation of deep learning input
Based on the results of both de novo search and MS-clust, NovoRank generates top two candidates.
The top two candidates are an initial point to train deep learning model.

A user can set the parameters in 'config_for_gen_top2.txt' file.
Parameter|Value|Explanation|Mandatory
---|---|---|---|
mgf_path|String|Path of a folder containing MS/MS spectra (MGF format).|Y|
denovo_result_csv|String|Path of the de novo search result CSV file (see Step 1. Preparation datasets).|Y|
db_result_csv|String|Path of the database search result CSV file (see Step 1. Preparation datasets).|N|
cluster_result_path|String|Path of the clustering result from MS-Cluster.|Y|
mgf_xcorr|String|Path of a folder containing MS/MS spectra for XCorr calculation (MGF format).|Y|
mgf_remove|String|Path of a folder containing MS/MS spectra to find internal fragment ions (MGF format).|Y|
precursor_search_ppm|Float|Precursor PPM tolerance.|Y|
elution_time|Integer|A total elution time in the mass spectrometry assay (minutes).|Y|
training|Boolean|If a user wants to train a model, set it True. Otherwise, set False (test only).|Y|
features_csv|String|Path of a result feature file as output.|Y|

Note that when training sets as "False", NovoRank ignores "db_result_csv".

```c
python gen_feature_top2_candidates.py config_for_gen_top2.txt
```

### Step 4. XCorr calculation
As a third-part, NovoRank uses XCorr value as an additional feature.

< Calculate XCorr using the following command line of CometX. >

```c
CometX.exe -X -Pcomet.params .\mgf_XCorr\*.mgf
``` 

### Step 5. The last step for training/test of NovoRank
Lastly, NovoRank takes three inputs: feature.csv and XCorr values obtained from Steps 3 and 4, respectively, as well as the MGF files.

A user can set the parameters in 'config_run_novorank.txt' file.
Parameter|Value|Explanation|Training or Test
---|---|---|---|
training|Boolean|If a user wants to train a model, set it True. Otherwise, set False (test only).|Both|
mgf_path|String|Path of a folder containing MS/MS spectra (MGF format).|Both|
mgf_xcorr|String|Path of the XCorr calculation TSV file.|Both|
features_csv|String|Path of the output of gen_feature_top2_candidates.py.|Both|
batch_size|Integer|Size of batch.|Both|
val_size|Float|The validation dataset ratio.|Training|
epoch|Integer|Size of epoch.|Training|
model_save_name|String|Save path and h5 file name for  trained model.|Training|
pre_trained_model|String|A path of pre-trained model h5 file.|Test|
result_name|String|Save path and CSV file name for test result.|Test|

"Both" means that it is used in both cases of training and Test.

```c
python run_novorank.py config_run_novorank.txt
```

#### Deep learning model for re-ranking.

The deep learning model only handles peptides with a maximum mass of 5000 Da and a length of 40 or less.

- Testing \
Using a pre-trained model, perform testing and output a single assigned peptide for each spectrum as the result.

- Training \
The deep learning model is trained based on the hyper-parameters set in the config_run_novorank.txt. \
The trained model is saved in the .h5 format as the output.
