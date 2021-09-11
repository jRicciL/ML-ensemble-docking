# Improving Structure-Based Virtual Screening with Ensemble Docking and Machine Learning

[![Hits](https://hits.seeyoufarm.com/api/count/incr/badge.svg?url=https%3A%2F%2Fgithub.com%2FjRicciL%2FML-ensemble-docking&count_bg=%2357C0C0&title_bg=%23555555&icon=jupyter.svg&icon_color=%23E7E7E7&title=hits&edge_flat=true)](https://hits.seeyoufarm.com)

## 🗂 About this repository

### 📙 **Content:**

- **Jupyter notebooks**, to follow the study's workflow,  
  required to reproduce the **results** and **figures** of the study. 
- **Helper Python and R scripts**.
- **Main datasets and complementary files**. 

### **📂 Protein directories:**

We evaluated **target-specific ML models** for *structure-based virtual screening*.  
The following **four proteins** were considered as case studies:

| # | Protein name | Directory  | UniProtKB  |
| - |- | - | ---- |
|1.  | **CDK2** | [`cdk2`](https://github.com/jRicciL/ML-ensemble-docking/tree/main/cdk2)| [P24941](https://www.uniprot.org/uniprot/P24941)  | 
| 2. | **FXa**  | [`fxa`](https://github.com/jRicciL/ML-ensemble-docking/tree/main/fxa)| [P00742](https://www.uniprot.org/uniprot/P00742) |  
| 3. | **EGFR** | [`egfr`](https://github.com/jRicciL/ML-ensemble-docking/tree/main/egfr)| [P00533](https://www.uniprot.org/uniprot/P00533) | 
| 4. | **HSP90** | [`hsp90`](https://github.com/jRicciL/ML-ensemble-docking/tree/main/hsp90)| [P07900](https://www.uniprot.org/uniprot/P07900) | 

Each **protein** directory (`cdk2`, `fxa`, `egfr`, `hsp90`) has the following structure:
1. `📂 1_Download_and_prepare_protein_ensembles`:  
    - *Download and prepare protein crystalographic structures from `PDB`*
2. `📂 2_Molecular_libraries`
    - *Download and prepare ligand molecules from benchmarking sets*
3. `📂 3_Protein_Ensembles_Analysis`
    - *Create and analyze the protein ensembles*
4. `📂 4_Ensemble_docking_results`
    - *Prepare and gather Ensemble Docking results*
5. `📂 5_Machine_Learning`
    - *Evaluate consensus strategies and ML classifiers through 30x4cv*
    
### Requirements:
#### 🐍 `Python`: Conda environment and required python libraries
- Install [conda](https://conda.io/projects/conda/en/latest/user-guide/install/index.html) or [miniconda](https://docs.conda.io/en/latest/miniconda.html).
- Create a conda environment from the `conda_environment.yml` file.  
```shell
conda env create -f conda_environment.yml
```

- The above will install all the python libraries used during our study.

#### 🔵 `R`: Libraries used
- Some of the analysis and plots were performed using `R` (version 4.0.3)
- The `R` libraries used here are listed at the top of each `R` script inside the `R_scripts` directory.


## 👥 Authors 

- **Joel Ricci-López**: *CICESE Research Center, Ensenada, México*
- **Sergio A. Aguila**: *CNyN, UNAM, Ensenada, México*
- **Michael K. Gilson**: *Skaggs School of Pharmacy and Pharmaceutical Sciences,  
UCSD, La Jolla, California, USA.*
- **Carlos A. Brizuela**: *CICESE Research Center, Ensenada, México*
  


## ✅ Acknowledgements

- LANCAD-UNAM-DGTIC-286 and PAPIIT-DGAPA-UNAM-IG200320grants
- **CAB**  and  **JRL**  acknowledge  the  support  of  CONACyT  under  grant  A1-S-20638
- **JRL**  was  supported  by  the  Programa  de  Doctorado  en  Nanociencias  at  CICESE  and  byCONACyT.
- Authors also thank to the **anonymous reviewers** for their comments and thoughtful suggestions, which substantially helped to improve the manuscript.

  
