# 2025_MMM_study

**Menstrual Cycle, Motivation and Mental Health Study (2025)**

Annalise Whines, Sara Mehrhof, Tory Eisenlohr-Moul & Camilla Nord

## Data

All data used in our study can be found in the `/data` directory. Participant IDs have been replaced with randomly assigned study IDs, and demographic information is limited to age and group.

### Processed data

The cleaned data resulting from our study can be found in the `/processed_dat` directory. Data processing includes applying pre-registered exclusion criteria to both questionnaire and task data. This directory also contains intermediate datasets created by the analysis scripts (`mhq_analysis.RDS` and `combined_params.RDS`), which are used in the change analyses.

### Model fits

Model fit objects produced by the model-based scripts in the `/analyses` directory, organised by cycle phase (`/menstrual` and `/ovulation`).

## Code

All code used to run analyses for our study can be found in the `/code` directory. Some scripts may call supporting functions and Stan scripts found in the `/functions` and `/stan` directories.

### Analyses

All code used for our main study analyses can be found in the `/analyses` directory. The analysis steps are structured as follows:

**1 Task – model based**
Model-based analyses of task data, run separately for the periovulatory and perimenstrual phases. Includes model fitting, checking for model convergence, model comparison and posterior predictive checks.

**2 Mental health analyses**
Comparison of mental health questionnaire measures across the periovulatory and perimenstrual phases.

**3 Task – main analyses**
Comparison of computational model parameters (effort sensitivity, reward sensitivity and acceptance bias) across the periovulatory and perimenstrual phases.

**4 Change analyses**
Associations between cycle-related changes in model parameters and changes in mental health measures.

The `/output` directory is empty and provides the folder structure for figures generated when running the scripts.
