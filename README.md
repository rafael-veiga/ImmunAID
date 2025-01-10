# ImmunAID
Code used in analyse for the paper:

**Patients with inflammation of unknown origin phenocopy the immune presentation of adult onset Still disease**

## Table of contents
* [Dataset description](#dataset-description)
* [Files description](#files-description)
* [Necessary repositories](#necessary-repositories)
## Dataset description
 disponible at Supplementary Spreadsheet 2.csv
 
 187 cases, 212 variables
 
| Colums Name |  type | description |
| :-----: | :------: | :------: |
| id | String | patient register |
| disease | categorical | presence of specific disease |
| age | integer | pacient age in years at sample colect |
| sex | categoric | sex of the patient |
| < "cell variant" > in < "cell population" > | Float number | estimated quantity of cells variant in an especific cell population |
  
## Files description
* R Files
  * **functions_pre.R** : set of functions used in all R scripts.
  * **Construct_dataset.R** : Script that integrate data sorces and create a dataset (not neccessary for the use of dataset provide).
  * **A1_pre_process_data_analise.R** : Script execute pre-process data (not neccessary for the use of dataset provide):
    * Remove missing
    * Transformation
    * Normalization
    * Imputation.
  * **A2_save_raw_data_csv.R** : Script execute transformations on provide dataset:
    * Transformation
    * Normalization
    * Imputation. 
  * **B3_Descriptive_table.R** : Construct Table 1:
  * **C3_Figs_create.R** : Construct all paper Figures and evaluate logistic regression immunologic marks odds ration. Split in 3 parts:
    * Base: load librarys, functions and definitions. (always need to run before the nexts parts).
    * Part 1: calculate magnitude of effect in each disiase for each immunological marker. the estimeted effect is avalieted by odds ration estimated by logistic regression.
    * Part 2: create all figures of the paper (it is necessary to run after Part 1 all python scripts before execute Part 2).
* Python Files
  * **analise4.py** : Execute estimation of diferent model generalization and behave.
  * **auc_curv.py** : Execute evaluation of models in relation to increase the number of relevant marks.
   
## Dependences
This code was tested on linux and windows machines. A standart computer is suficient to run all analisys and should not take more than hours for each step.
The follow softwares and packadges are necessary:
* **R** : version 4.4.0
* **python** : version 3.9.18
 * **numpy** : version 1.23.5
 * **pandas** : version  1.5.3 
 * **scikit-learn** : version 1.2.1

## Instruction To Run
* **1 - Download data** : from the repository http://flowrepository.org/id/FR-FCM-Z662 download the file **data_raw.rds** and put inside **pos_data** folder 
* **2 - Run A1_pre_process_data_analise.R** : execute all file content in R.
* **3 - Run A2_save_raw_data_csv.R** : execute all file content in R.
* **4 - Run B3_Descriptive_table.R** : execute all file content in R.
* **5 - Run C3_Figs_create.R** : execute in R **Base** and **Part1**.
* **6 - run analise4.py** : execupe in python all content
* **7 - Run C3_Figs_create.R** : execute in R **Base** and **Part2**.
