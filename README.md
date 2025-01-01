# NovoRank: Machine Learning-Based Post-Processing Tool for *De Novo* Peptide Sequencing

NovoRank is a post-processing tool designed to improve the accuracy of *de novo* peptide sequencing in proteomics. Unlike database-dependent methods, *de novo* sequencing derives peptide sequences directly from tandem mass spectrometry (MS/MS) data, enabling the discovery of novel peptides. However, reliance on incomplete scoring functions often leads to incorrect identifications. NovoRank addresses this by re-ranking candidate peptides to recover correct identifications, enhancing both precision and recall. It is compatible with any *de novo* sequencing software, and by reassigning the optimal peptide, NovoRank offers a robust solution for overcoming the noise and ambiguities inherent in MS/MS data.

For detailed insights behind NovoRank, refer to the [**NovoRank paper**](https://pubs.acs.org/doi/10.1021/acs.jproteome.4c00300).

<br>

<p align="center">
    <img src="assets/img.png" width="700"\>
</p>
<hr>

## Overview
#### Code structure
``` Unicode
NovoRank
  │ 
  ├─── generate_candidates_and_extract_features.py
  │
  ├─── run_novorank.py
  │ 
  ├─── src
  │     │     
  │     ├── features
  │     │     └── featureprocessor.py: Functions for feature calculation
  │     │       
  │     ├── loader
  │     │     └── dataloader.py: Functions for data loading
  │     │
  │     ├── model
  │     │     ├── base_model.py: Base model structure
  │     │     ├── inference.py: Functions for model inference
  │     │     ├── preprocess.py: Functions for data preprocessing in modeling
  │     │     └── train.py: Functions for model training
  │     │ 
  │     └─── utils
  │           ├── config_first.py: Command-line argument parsing and configuration file loading for generate_candidates_and_extract_features.py
  │           ├── config_second.py: Command-line argument parsing and configuration file loading for run_novorank.py
  │           ├── process.py: Functions for performing data processing
  │           └── utils.py: Utility functions providing support
  │
  ├─── models: Trained models (may include models saved at each epoch)
  │
  ├─── pretrained: Pretrained NovoRank models (Casanovo, PEAKS, pNovo3) in .h5 format
  │
  └─── software
        ├── CometX: XCorr calculation software (in-house software)
        └── MSCluster: Spectral clustering software
```

## Datasets
All datasets used in this work are available for download from [**Zenodo**](https://zenodo.org/records/14046459).  

To use NovoRank, a user MUST refer to the `README.md` in the [./data](https://github.com/jangho721/NovoRank/tree/main/data) directory, where sample data has also been provided.

## Configuration
  
  The [config.yaml](https://github.com/jangho721/NovoRank/blob/main/config.yaml) is used to set up the parameters and initial configurations required to run NovoRank. It contains default values, and descriptions for each option are provided as comments within it.

## Requirements
⦁ To install the required Python packages:
1. Clone the repository or download the code.
2. Create and activate an Anaconda virtual environment:
```c
conda create -n [NAME] python==3.9
conda activate [NAME]
```
3. Install the dependencies listed in `requirements.txt`:
```c
pip install -r requirements.txt
```
> **Note:**  
> NovoRank was implemented using **Python 3.9** and utilizes the [**DeepLC**](https://github.com/compomics/DeepLC) package, which is included in the `requirements.txt`.
<br>

⦁ Software
- [**MS-Cluster**](http://proteomics.ucsd.edu/software-tools/ms-clusterarchives) ([Download](https://github.com/jangho721/NovoRank/tree/main/software/MSCluster))
- [**CometX**](https://github.com/jangho721/NovoRank/tree/main/software/CometX) (In-house software modified to calculate XCorr, based on Comet software)

## How to Use
For the description of the datasets required to execute Steps 2 and 4, refer to [Essential Data for Using NovoRank](https://github.com/jangho721/NovoRank/blob/main/data/README.md#essential-data-for-using-novorank).  

- Step 1. Spectral clustering using MS-Cluster
```c
MSCluster.exe --list [PATH] --output-name CLUSTERS --mixture-prob 0.01 --fragment-tolerance 0.02 --assign-charges
```
> **Note:**  
> For detailed instructions on using MS-cluster, refer to the [manual](https://github.com/jangho721/NovoRank/blob/main/software/MSCluster/manual.pdf).
<br>

- Step 2. Generate two candidates and extract features

The parameter that controls the training process or the inference process is specified in the [config.yaml](https://github.com/jangho721/NovoRank/blob/main/config.yaml).  
```c
python generate_candidates_and_extract_features.py --search_ppm [PRECURSOR_TOLERANCE] --elution_time [ELUTION_TIME_MIN]
```
> **Note:**  
> To check the available options and their descriptions, run the command `python generate_candidates_and_extract_features.py -h`.
<br>

- Step 3. Xcorr calculation using CometX  
```c
CometX.exe -X -Pcomet.params [PATH]\*.mgf
```
> **Note:**  
> [PATH] is the directory containing MGF files for xcorr calculation.  
> `*_xcorr.tsv` files will be generated.
<br>

- Step 4. Training & Inference of the NovoRank  
```c
python run_novorank.py
```
> **Note:**  
> To check the available options and their descriptions, run the command `python run_novorank.py -h`.
<br>

1. **Training**  
To use NovoRank, users **must train a model tailored to their dataset** - it is **recommended** to use a customized model based on the *de novo* search software used.
The trained model is saved in the `./models/` directory in `.h5` format. Additionally, checkpoint models trained at each epoch can be saved.
<br>

2. **Inference**  
Inference can be performed using the pre-trained model. The pre-trained models for testing, created using three types of *de novo* search software (Casanovo, PEAKS, pNovo3), are located in the `./pretrained/` directory.
<br>
The deep learning model only handles peptides with a maximum mass of 5000 Da and a length of 40 or less.

## Results
<pre>
- The results_top1.csv file is generated at the ./data/interim location.
  (The save location and result file name can be changed in the config.yaml file)
- The NovoRank results are to output a single assigned peptide for each spectrum.
</pre>

## Credits
NovoRank is created by <a href="https://jangho721.github.io/" target="_blank">Jangho Seo</a>, Seunghyuk Choi, and Eunok Paek at the Hanyang University.

## Citation
```bibTeX
@article{sep2024novorank,
  title = {NovoRank: Refinement for De Novo Peptide Sequencing Based on Spectral Clustering and Deep Learning},
  shorttitle = {NovoRank},
  author = {Seo, Jangho and Choi, Seunghyuk and Paek, Eunok},
  journal={Journal of Proteome Research},
  year={2024}
}
```

## License
<pre>
- NovoRank © 2024 is licensed under Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International.
  This license requires that reusers give credit to the creator. It allows reusers to distribute, 
  remix, adapt, and build upon the material in any medium or format, for noncommercial purposes only. 
  If others modify or adapt the material, they must license the modified material under identical terms.
</pre>

## Contact
If you have any questions, feel free to [open an issue](https://github.com/jangho721/NovoRank/issues/new) or contact [Jangho Seo](https://jangho721.github.io/) or any of the contributors listed above.
