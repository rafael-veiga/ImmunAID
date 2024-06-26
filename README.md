# ImmunAID
Code used in analyse for the paper:

**Patients with inflammation of unknown origin phenocopy the immune presentation of adult onset Still disease**

## Table of contents
* [Dataset description](#dataset-description)
* [Files description](#files-description)
* [Necessary repositories](#necessary-repositories)
## Dataset description
 disponible at /pos_data/ImmunAID.csv
 
 310 cases, 505 variables
 
| Colums Name |  type | description |
| :-----: | :------: | :------: |
| id | integer | patient register |
| disease | categorical | presence of specific disease |
| age | integer | pacient age in years at sample colect |
| sex | categoric | sex of the patient |
| fil_dis | boolean | false if the patient has a disease not belong to this study.
| < "cell variant" > in < "cell population" > | Float number | estimated quantity of cells variant in an especific cell population |
  
## Files description
* R Files
  * **functions_pre.R** : set of functions used in all R scripts.
  * **Construct_dataset.R** : Script that integrate data sorces and create a dataset (not neccessary for the use of dataset provide).
  * **A2_pre_process_data_analise.R** : Script execute pre-process data:
    * Remove missing
    * Transformation
    * Normalization
    * Imputation.
  * **B3_Descriptive_table.R** : Construct Table 1:
  * **C3_Figs_create.R** : Construct all paper Figures and evaluate logistic regression immunologic marks odds ration. Divide in 3 parts:
    * Base: load librarys, functions and definitions. (always need to run before the nexts parts).
    * Part 1: calculate magnitude of effect in each disiase for each immunological marker. the estimeted effect is avalieted by odds ration estimated by logistic regression.
    * Part 2: create all figures of the paper (it is necessary to run after Part 1 all python scripts before execute Part 2).
* Python Files
  * **analise4.py** : Execute estimation of diferent model generalization and behave.
  * **auc_curv.py** : Execute evaluation of models in relation to increase the number of relevant marks.
   
## Necessary repositories
* **R** : version 4.4.0
*  **python** : version 3.9.18
*  **numpy** : version 1.23.5
* **pandas** : version  1.5.3 
* **scikit-learn** : version 1.2.1
