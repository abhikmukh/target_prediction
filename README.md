# Python repositiory for predicting drug targets from SMILES string

## Overview
This is a python repository to predict probable target for a small molecule. This work is based on work from these two papers "prediction_of_bayessian_target_on_multiple_database" by Jenkins et.al (J. Chem. Inf. Model. 2006, 46, 1124-1133) and "Searching Chemical Space with the Bayesian Idea Generator" (J. Chem. Inf. Model. 2009, 49, 2211–2220) by van Hoorn et.al. 
In this work, I have used a subset of target-ligand dataset from BindingDB. This repo generates multi category Naive Bayesian models based from the dataset and load the models in a locally hosted PostGreSQL server. A smiple web application has been created to run the python script. The web application returns top 5 predicted uniprot target ids for an input ligand. 

## Requirements
To run this repository create a conda environemnt with this command conda env create --name envname --file=environments.yml


## How to run this project
Once the conda enviroment is created run the process.py inside the conda environment. It will create the model and load the model in a local PostgreSQL database. Next run the app.py to start the web server. Once the web server is started please follwo the instructions in the web page. 


![web screenshot](/screenshots/webpage_screeshot.PNG?raw=true)

## Database schema
![database schema](/screenshots/schema.PNG?raw=true)

## License information
Apache License 2.0