# FXa jupyter notebooks

## **Directory structure:**

### 1. **`Download_and_prepare_protein_ensembles`:**

- Download, model, and prepare <mark style='background-color: #FFF2CD'>crystallographic protein</mark> structures from the Protein Data Bank.

    
| # | Description | Notebook  | View |
| - |- | - | ---- |
| 1 | Download `PDB` structures | 📙 [`1_Download_crystal_structures_from_PDB.ipynb`](https://github.com/jRicciL/ML-ensemble-docking/tree/main/fxa/1_Download_and_prepare_protein_ensembles/1_Download_crystal_structures_from_PDB.ipynb) | [![View the notebook](https://img.shields.io/badge/render-nbviewer-orange.svg)](https://nbviewer.jupyter.org/github/jRicciL/ML-ensemble-docking/blob/main/fxa/1_Download_and_prepare_protein_ensembles/1_Download_crystal_structures_from_PDB.ipynb) |
| 2 | Get `PDB` metadata | 📙 [`2_Get_PDB_structures_metadata.ipynb`](https://github.com/jRicciL/ML-ensemble-docking/tree/main/fxa/1_Download_and_prepare_protein_ensembles/2_Get_PDB_structures_metadata.ipynb) | [![View the notebook](https://img.shields.io/badge/render-nbviewer-orange.svg)](https://nbviewer.jupyter.org/github/jRicciL/ML-ensemble-docking/blob/main/fxa/1_Download_and_prepare_protein_ensembles/2_Get_PDB_structures_metadata.ipynb) |
| 3 | Model missing loops | 📙 [`3_Model_structures_using_Modeller.ipynb`](https://github.com/jRicciL/ML-ensemble-docking/tree/main/fxa/1_Download_and_prepare_protein_ensembles/3_Model_structures_using_Modeller.ipynb) | [![View the notebook](https://img.shields.io/badge/render-nbviewer-orange.svg)](https://nbviewer.jupyter.org/github/jRicciL/ML-ensemble-docking/blob/main/fxa/1_Download_and_prepare_protein_ensembles/3_Model_structures_using_Modeller.ipynb) |
| 4 | Prepare protein conformations | 📙 [`4_Prepare_proteins_using_pdb4amber.ipynb`](https://github.com/jRicciL/ML-ensemble-docking/tree/main/fxa/1_Download_and_prepare_protein_ensembles/4_Prepare_proteins_using_pdb4amber.ipynb) | [![View the notebook](https://img.shields.io/badge/render-nbviewer-orange.svg)](https://nbviewer.jupyter.org/github/jRicciL/ML-ensemble-docking/blob/main/fxa/1_Download_and_prepare_protein_ensembles/4_Prepare_proteins_using_pdb4amber.ipynb) |
| 5 | Get `COCRYS` molecules | 📙 [`5_Get_cocrystalized_molecules_from_PDB.ipynb`](https://github.com/jRicciL/ML-ensemble-docking/tree/main/fxa/1_Download_and_prepare_protein_ensembles/5_Get_cocrystalized_molecules_from_PDB.ipynb) | [![View the notebook](https://img.shields.io/badge/render-nbviewer-orange.svg)](https://nbviewer.jupyter.org/github/jRicciL/ML-ensemble-docking/blob/main/fxa/1_Download_and_prepare_protein_ensembles/5_Get_cocrystalized_molecules_from_PDB.ipynb) |



### 2. **`Molecular_libraries`:**

- Download and prepare the <mark style='background-color: #FFF2CD'>molecular libraries</mark>.
    
| # | Description | Notebook  | View |
| - |- | - | ---- |
| 1 | Download molecular libraries |📙 [`1_Get_Molecular_libraries.ipynb`](https://github.com/jRicciL/ML-ensemble-docking/tree/main/fxa/2_Molecular_libraries/1_Get_Molecular_libraries.ipynb) | [![View the notebook](https://img.shields.io/badge/render-nbviewer-orange.svg)](https://nbviewer.jupyter.org/github/jRicciL/ML-ensemble-docking/blob/main/fxa/2_Molecular_libraries/1_Get_Molecular_libraries.ipynb) |
| 2 | Load molecules with `rdkit` |📙 [`2_Loading_molecules_from_db_with_rdkit.ipynb`](https://github.com/jRicciL/ML-ensemble-docking/tree/main/fxa/2_Molecular_libraries/2_Loading_molecules_from_db_with_rdkit.ipynb) | [![View the notebook](https://img.shields.io/badge/render-nbviewer-orange.svg)](https://nbviewer.jupyter.org/github/jRicciL/ML-ensemble-docking/blob/main/fxa/2_Molecular_libraries/2_Loading_molecules_from_db_with_rdkit.ipynb) |
| 3 | Find duplicates |📙 [`3_Comparing_Molecules_among_Molecular_Libraries.ipynb`](https://github.com/jRicciL/ML-ensemble-docking/tree/main/fxa/2_Molecular_libraries/3_Comparing_Molecules_among_Molecular_Libraries.ipynb) | [![View the notebook](https://img.shields.io/badge/render-nbviewer-orange.svg)](https://nbviewer.jupyter.org/github/jRicciL/ML-ensemble-docking/blob/main/fxa/2_Molecular_libraries/3_Comparing_Molecules_among_Molecular_Libraries.ipynb) |
| 4 |EDA of molecular libraries |📙 [`4_Data_Visualization_of_molecular_libraries.ipynb`](https://github.com/jRicciL/ML-ensemble-docking/tree/main/fxa/2_Molecular_libraries/4_Data_Visualization_of_molecular_libraries.ipynb) | [![View the notebook](https://img.shields.io/badge/render-nbviewer-orange.svg)](https://nbviewer.jupyter.org/github/jRicciL/ML-ensemble-docking/blob/main/fxa/2_Molecular_libraries/4_Data_Visualization_of_molecular_libraries.ipynb) |

### 3. **`Protein_Ensembles_Analysis`:**

- Prepare and analyze protein conformational ensembles.
    
| # | Description | Notebook  | View |
| - |- | - | ---- |
| 1 | Analize protein ensembles | 📙 [`1_Structural_analysis_of_Protein_Ensembles.ipynb`](https://github.com/jRicciL/ML-ensemble-docking/tree/main/fxa/3_Protein_Ensembles_Analysis/1_Structural_analysis_of_Protein_Ensembles.ipynb) | [![View the notebook](https://img.shields.io/badge/render-nbviewer-orange.svg)](https://nbviewer.jupyter.org/github/jRicciL/ML-ensemble-docking/blob/main/fxa/3_Protein_Ensembles_Analysis/1_Structural_analysis_of_Protein_Ensembles.ipynb) |

### 4. **`Ensemble_docking_results`:**

- Analyze ensemble docking results.
    
| # | Description | Notebook  | View |
| - |- | - | ---- |
| 1 | `SMINA/Vinardo` docking params. | 📙 [`1_Docking_Parameters_used.ipynb`](https://github.com/jRicciL/ML-ensemble-docking/tree/main/fxa/4_Ensemble_docking_results/1_Docking_Parameters_used.ipynb) | [![View the notebook](https://img.shields.io/badge/render-nbviewer-orange.svg)](https://nbviewer.jupyter.org/github/jRicciL/ML-ensemble-docking/blob/main/fxa/4_Ensemble_docking_results/1_Docking_Parameters_used.ipynb) |
| 2 | Load Enseble Docking resutls | 📙 [`2_Ensemble_Docking_Results.ipynb`](https://github.com/jRicciL/ML-ensemble-docking/tree/main/fxa/4_Ensemble_docking_results/2_Ensemble_Docking_Results.ipynb) | [![View the notebook](https://img.shields.io/badge/render-nbviewer-orange.svg)](https://nbviewer.jupyter.org/github/jRicciL/ML-ensemble-docking/blob/main/fxa/4_Ensemble_docking_results/2_Ensemble_Docking_Results.ipynb) |
| 3 | Docking results per conformation| 📙 [`3_Ensemble_Docking_results_per_conformation.ipynb`](https://github.com/jRicciL/ML-ensemble-docking/tree/main/fxa/4_Ensemble_docking_results/3_Ensemble_Docking_results_per_conformation.ipynb) | [![View the notebook](https://img.shields.io/badge/render-nbviewer-orange.svg)](https://nbviewer.jupyter.org/github/jRicciL/ML-ensemble-docking/blob/main/fxa/4_Ensemble_docking_results/3_Ensemble_Docking_results_per_conformation.ipynb) |

### 5. **`Machine_Learning`:**

- Apply and compare `Consensus strategies` and `Machine learning` classifiers through 30 repetitions of 4-fold cross-validation ($30 \times 4 cv$), conformational selection strategies, and `Recursive Feature Elimination`.
    
| # | Description | Notebook  | View |
| - |- | - | ---- |
| 1 | Hyperparameter tuning | 📙 [`1_Hyperparameter_tuning_stage.ipynb`](https://github.com/jRicciL/ML-ensemble-docking/tree/main/fxa/5_Machine_Learning/1_Hyperparameter_tuning_stage.ipynb) | [![View the notebook](https://img.shields.io/badge/render-nbviewer-orange.svg)](https://nbviewer.jupyter.org/github/jRicciL/ML-ensemble-docking/blob/main/fxa/5_Machine_Learning/1_Hyperparameter_tuning_stage.ipynb) |
| 2a | __*30x4cv*__: All confs., tuned Hyprms. | 📙 [`2_30x4CV_analysis.ipynb`](https://github.com/jRicciL/ML-ensemble-docking/tree/main/fxa/5_Machine_Learning/2_30x4CV_analysis.ipynb) | [![View the notebook](https://img.shields.io/badge/render-nbviewer-orange.svg)](https://nbviewer.jupyter.org/github/jRicciL/ML-ensemble-docking/blob/main/fxa/5_Machine_Learning/2_30x4CV_analysis.ipynb) |
| 2b | __*30x4cv*__: All confs., default Hyprms. |📙 [`2_30x4CV_analysis-DEFAULT-HYPERPARAMETERS.ipynb`](https://github.com/jRicciL/ML-ensemble-docking/tree/main/fxa/5_Machine_Learning/2_30x4CV_analysis-DEFAULT-HYPERPARAMETERS.ipynb) | [![View the notebook](https://img.shields.io/badge/render-nbviewer-orange.svg)](https://nbviewer.jupyter.org/github/jRicciL/ML-ensemble-docking/blob/main/fxa/5_Machine_Learning/2_30x4CV_analysis-DEFAULT-HYPERPARAMETERS.ipynb) |
| 2c | __*30x4cv*__: All confs., Ligand Efficiency scores |📙 [`2_30x4CV_analysis-LIG-EFF.ipynb`](https://github.com/jRicciL/ML-ensemble-docking/tree/main/fxa/5_Machine_Learning/2_30x4CV_analysis-LIG-EFF.ipynb) | [![View the notebook](https://img.shields.io/badge/render-nbviewer-orange.svg)](https://nbviewer.jupyter.org/github/jRicciL/ML-ensemble-docking/blob/main/fxa/5_Machine_Learning/2_30x4CV_analysis-LIG-EFF.ipynb) |
| 2d | Repeated CV versions| 📙 [`30x4CV_vs_NestedCV.ipynb`](https://github.com/jRicciL/ML-ensemble-docking/tree/main/fxa/5_Machine_Learning/30x4CV_vs_NestedCV.ipynb) | [![View the notebook](https://img.shields.io/badge/render-nbviewer-orange.svg)](https://nbviewer.jupyter.org/github/jRicciL/ML-ensemble-docking/blob/main/fxa/5_Machine_Learning/30x4CV_vs_NestedCV.ipynb) |
| 3 | *y-scrambling* validation | 📙 [`3_y_scambling_analysis.ipynb`](https://github.com/jRicciL/ML-ensemble-docking/tree/main/fxa/5_Machine_Learning/3_y_scambling_analysis.ipynb) | [![View the notebook](https://img.shields.io/badge/render-nbviewer-orange.svg)](https://nbviewer.jupyter.org/github/jRicciL/ML-ensemble-docking/blob/main/fxa/5_Machine_Learning/3_y_scambling_analysis.ipynb) |
| 4 | Feature ranking with **RFE** | 📙 [`4_Recursive_Feature_Elimination.ipynb`](https://github.com/jRicciL/ML-ensemble-docking/tree/main/fxa/5_Machine_Learning/4_Recursive_Feature_Elimination.ipynb) | [![View the notebook](https://img.shields.io/badge/render-nbviewer-orange.svg)](https://nbviewer.jupyter.org/github/jRicciL/ML-ensemble-docking/blob/main/fxa/5_Machine_Learning/4_Recursive_Feature_Elimination.ipynb) |
| 5 | Conf. selection using *k* confs.| 📙 [`5_Conformational_selection_using_k_confs.ipynb`](https://github.com/jRicciL/ML-ensemble-docking/tree/main/fxa/5_Machine_Learning/5_Conformational_selection_using_k_confs.ipynb) | [![View the notebook](https://img.shields.io/badge/render-nbviewer-orange.svg)](https://nbviewer.jupyter.org/github/jRicciL/ML-ensemble-docking/blob/main/fxa/5_Machine_Learning/5_Conformational_selection_using_k_confs.ipynb) |
| 6 | __*30x4cv*__ usign *k* confs. |📙 [`6_Conformational_selection_Plots.ipynb`](https://github.com/jRicciL/ML-ensemble-docking/tree/main/fxa/5_Machine_Learning/6_Conformational_selection_Plots.ipynb) | [![View the notebook](https://img.shields.io/badge/render-nbviewer-orange.svg)](https://nbviewer.jupyter.org/github/jRicciL/ML-ensemble-docking/blob/main/fxa/5_Machine_Learning/6_Conformational_selection_Plots.ipynb) |
| 7 | cMDS plots | 📙 [`7_Best_Confs_MDS.ipynb`](https://github.com/jRicciL/ML-ensemble-docking/tree/main/fxa/5_Machine_Learning/7_Best_Confs_MDS.ipynb) | [![View the notebook](https://img.shields.io/badge/render-nbviewer-orange.svg)](https://nbviewer.jupyter.org/github/jRicciL/ML-ensemble-docking/blob/main/fxa/5_Machine_Learning/7_Best_Confs_MDS.ipynb) |
| 8 | Top ranked molecules per method | 📙 [`8_Top_ranking_molecules_by_method.ipynb`](https://github.com/jRicciL/ML-ensemble-docking/tree/main/fxa/5_Machine_Learning/8_Top_ranking_molecules_by_method.ipynb) | [![View the notebook](https://img.shields.io/badge/render-nbviewer-orange.svg)](https://nbviewer.jupyter.org/github/jRicciL/ML-ensemble-docking/blob/main/fxa/5_Machine_Learning/8_Top_ranking_molecules_by_method.ipynb) |