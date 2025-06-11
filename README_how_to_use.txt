Data and code to replicate the analyses reported in the paper:
Ez-zizi, A., Farrell, S., Leslie, D., Malhotra, G., and Ludwig, J.H.C. Reinforcement learning under uncertainty: expected versus unexpected uncertainty and state versus reward uncertainty.

All scripts were written by Dr Adnane Ez-zizi (a.ez-zizi@uos.ac.uk)

---------
-- Data
---------

The datasets generated for each participant in each experiement (1 or 2) and condition (control versus test in Exp 1 
and state versus reward in Exp 2) can be found in the "Data" folder. The folder also contains datasets in both 
csv and rds formats, which record data from all participants in each experiment. These datasets are named "Data_all".    

---------
-- Code
---------

The package is made up of folders containing scripts for running the behavioural experiments, fitting the 
computationnal models to the data and simulating them, and for reproducing the statisical analyses reported in the 
paper. 

I) Experiment_scripts 
--------------------

The scripts that were used to run both experiments are stored in the folder "Experiment_scripts". Please read the 
README file in the folder Experiment2, which explains in detail how to run the second experiment. The intruction sheets
are also included. 


II) Computational_modelling
---------------------------

The scripts for fitting the different computational models are provided in the folder "Computational_modelling". The 
folder should be easy to navigate as it is divided by experiment and condition. the folder ifit is the package used to 
run simulated annealing (downloaded from http://ifit.mccode.org/). The two files that you will need to run for each model
are:

- Fitting_file.m: This will fit the model to the data and estimate the model parameters (e.g. learning rate, exploration).
  After running the script, you will get the results both as .mat and .txt files (FittedParameters.txt and Fit<CONDITION>.mat).
  Please note that the models can take several hours to fit especially the BISAW model. You could lower the parameter Ns 
  and/or the number of random parameter initiations to speed up the runs, but the fit precision will be lower.

- SimulationFile_fit.m: This will simulated the fitted model on the same observations that the fitted participant encountered.
  This is what we used for example to make Figure 10. 

III) Statistical_analyses
------------------------

The folder "Statistical_analyses" contains R scripts for reproducing the results of the GAMM analyses (e.g. Tables 1, 2 and 5; 
Figures 6, 9 and 14). The scripts start with a short description that explains what results they reporoduce. For Experiment 2, 
please follow the order suggested by the prefix numbers when executing the sctipts. 
  


