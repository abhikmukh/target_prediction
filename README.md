# Python repositiory for predicting drug targets from SMILES string

## Overview
This is a python repository to predict targets from SMILES. This work is based "prediction_of_bayessian_target_on_multiple_database" by Jenkins et.al (J. Chem. Inf. Model. 2006, 46, 1124-1133) and "Searching Chemical Space with the Bayesian Idea Generator" (J. Chem. Inf. Model. 2009, 49, 2211–2220) by van Hoorn et.al. 
In this work, I have used a subset of target-ligand dataset from BindingDB. This repo generates multi category Naive Bayesian models based from the dataset and load the models in a locally hosted PostGreSQL server. A smiple web application has been created to run the python script. The web application returns top 5 predicted uniprot target ids for a SMILES of a ligand. 


