beg = clock; % starting time

% Add the path to the ifit package folder (MODIFY AS NEEDED)
addpath(genpath('./Computational_modelling/ifit'));

paramVector = zeros(5,1);
paramVector(1) = 200; % N - total nb of steps in the simulation
paramVector(2) = 0.7; % rho - level of confidence (in perception)
paramVector(3) = 0.02; % alpha - learning rate
paramVector(4) = 0.05; % T - temperature
paramVector(5) = 10; % r0 - reward value for correct responses

[N,rho,alpha0,T0,r0] = DispParams(paramVector);
J = 31; % number of subjects
numOfParams = 3; % number of parameters to estimate
Ns = 1000; % nb of simulations to use to estimate the action likelihoods - Ns (1000)
numOfParamInit = 30; % Num of repetitions we run the fitting procedure for each model-subject - we randomly generated 30 initial parameter values
NumOfWorkers = 16; % Num of cores to use in the foreach loop

%create parallel pool of workers on the local node
%Ensure that this is the same number as what you requested from the scheduler
pool = parpool('local', NumOfWorkers);

% Define filenames of result file (record results for each parameter value averaged over the number of simulation runs):
datafilename1 = strcat('./Results/FittedParameters.txt'); % name of data file to write to
datafilepointer1 = fopen(datafilename1,'wt'); % open ASCII file for writing
fprintf(datafilepointer1,'Subj rho_fit alpha_fit T_fit LogML AIC BIC\r\n'); % Heading of the data file

%%% Initial parameters for the fitting procedure
rng('default');
rng(1); % set seed to be able to replicate the results
alpha0_vals = rand(numOfParamInit,1);
T0_vals = rand(numOfParamInit,1);
rho0_vals=(1+rand(numOfParamInit,1))/2; % to make sure rho is initialised between 0.5 and 1 

%%% Matrix containing fitted_alpha, fitted_T and log(ML) values for each starting value of alpha and T for one subject
N_alpha0 = length(alpha0_vals); % number of starting values for alpha
N_T0 = length(T0_vals); % number of starting values for T
N_rho0 = length(rho0_vals); % number of starting values for rho
Matj = zeros(N_alpha0,7);

% Matrix containing the best fitted values with the corresponding logML for each subject
MatFittedParams = zeros(J,7);
MatFittedParams(:,1) = (1:J)';  

beg_subj = clock;

for subj = 1:J
    
    % Define filenames of result file (record results for each parameter value averaged over the number of simulation runs):
    datafilename2 = strcat('./Results/FittedParameters_subj',num2str(subj+2),'.txt'); % name of data file to write to
    datafilepointer2 = fopen(datafilename2,'wt'); % open ASCII file for writing
    fprintf(datafilepointer2,'rho0 alpha0 T0 LogML rho_fit alpha_fit T_fit\r\n'); % Heading of the data file 
    
    % Change to the path containing the behavioural data files
    Dataj = ['./Data/StateUncertaintyPhase2_',num2str(subj+2),'.txt'];
    
    parfor ii = 1:size(alpha0_vals,1)

        tic; % tic toc to measure run time
        
        alpha0 = alpha0_vals(ii);
        T0 = T0_vals(ii);
        rho0 = rho0_vals(ii);

        [alpha_fit,T_fit,rho_fit,LogML] = PWRL_fitting([alpha0,T0,rho0],r0,N,Ns,Dataj);
        
        % store the fitted values
        Matj(ii,:) = [rho0,alpha0,T0,rho_fit,alpha_fit,T_fit,LogML];
        
        toc  
        
    end
    
    % Record Matj into a text file
    for nn = 1:size(alpha0_vals,1)
        
        % Write trial result to file:   
        fprintf(datafilepointer2,'%1.3f %1.3f %1.3f %1.4f %1.3f %1.3f %1.3f\r\n', ...
                Matj(nn,1), ...
                Matj(nn,2), ...
                Matj(nn,3), ...
                Matj(nn,7), ...
                Matj(nn,4), ...
                Matj(nn,5), ...
                Matj(nn,6));

    end
    
    % get the set (rho_fit,alpha_fit,T_fit) that has the highest logML
    [logML_max, rowsOfMaxes] = max(Matj(:,7));
    
    rho_bestfit = Matj(rowsOfMaxes,4); % best fitted value for T
    alpha_bestfit = Matj(rowsOfMaxes,5); % best fitted value for alpha
    T_bestfit = Matj(rowsOfMaxes,6); % best fitted value for T
    [aic, bic] = aicbic(logML_max, numOfParams, N); % AIC and BIC scores
    
    MatFittedParams(subj,2) = rho_bestfit;
    MatFittedParams(subj,3) = alpha_bestfit;
    MatFittedParams(subj,4) = T_bestfit;
    MatFittedParams(subj,5) = logML_max;
    MatFittedParams(subj,6) = aic;
    MatFittedParams(subj,7) = bic;
    
    % Write trial result to file:   
    fprintf(datafilepointer1,'%i %1.6f %1.6f %1.6f %1.4f %1.2f %1.2f\r\n', ...
            subj+2, ...
            rho_bestfit, ...
            alpha_bestfit, ...
            T_bestfit, ...
            logML_max, ...
            aic, ...
            bic); 
        
    myfilename = ['./Results/FittingFile_subj_',num2str(subj+2),'.mat'];
    save(myfilename);
    
    % Output message to have a feel of the running progress
    ending_subj = clock;   
    fprintf('****Subject: %d - Elapsed time: %1.2f\n',[subj,etime(ending_subj,beg_subj)]);

end

% save fitting results in a '.mat' file
myfilename = './Results/FitTest.mat'; 
save(myfilename,'MatFittedParams');

% Wrapping up
delete(pool)
fclose('all');
ending = clock;
total_duration = etime(ending,beg);



