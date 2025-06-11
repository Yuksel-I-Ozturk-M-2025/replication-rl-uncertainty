beg = clock;

% Add the path to the ifit package folder (MODIFY AS NEEDED)
addpath(genpath('./Computational_modelling/ifit'));

paramVector = zeros(5,1);
paramVector(1) = 200; % N - total nb of steps in the simulation
paramVector(2) = 0.65; % rho - level of confidence (in perception)
paramVector(3) = 0.02; % alpha - learning rate
paramVector(4) = 0.05; % T - temperature
paramVector(5) = 2.5; % r0 - reward value for correct responses

[N,rho,alpha0,T0,r0] = DispParams(paramVector);
J = 12; % number of subjects
numOfParams = 2; % number of parameters to estimate
Ns = 1000; % nb of simulations to use to estimate the action likelihoods - Ns (1000)
numOfParamInit = 30; % Num of repetitions we run the fitting procedure for each model-subject - we randomly generated 30 initial parameter values
NumOfWorkers = 16; % Num of cores to use in the foreach loop

%create parallel pool of workers on the local node
%Ensure that this is the same number as what you requested from the scheduler 
pool = parpool('local', NumOfWorkers);

% Define filenames of result file (record results for each parameter value averaged over the number of simulation runs):
datafilename1 = strcat('./Results/FittedParameters.txt'); % name of data file to write to
datafilepointer1 = fopen(datafilename1,'wt'); % open ASCII file for writing
fprintf(datafilepointer1,'Subj alpha_fit T_fit LogML AIC BIC\r\n'); % Heading of the data file 

%%% Initial parameters for the fitting procedure
rng('default');
rng(1); % set seed to be able to replicate the results
alpha0_vals = rand(numOfParamInit,1);
T0_vals = rand(numOfParamInit,1);

%%% Matrix containing fitted_alpha, fitted_T and log(ML) values for each starting value of alpha and T for one subject
N_alpha0 = length(alpha0_vals); % number of starting values for alpha
N_T0 = length(T0_vals); % number of starting values for T
Matj = zeros(N_alpha0,5);

% Matrix containing the best fitted values with the corresponding logML for each subject
MatFittedParams = zeros(J,6);
MatFittedParams(:,1) = (1:J)'; 

beg_subj = clock;

for subj = 1:J
    
    % Define filenames of result file (record results for each parameter value averaged over the number of simulation runs):
    datafilename2 = strcat('./Results/FittedParameters_subj',num2str(subj+12),'.txt'); % name of data file to write to
    datafilepointer2 = fopen(datafilename2,'wt'); % open ASCII file for writing
    fprintf(datafilepointer2,'alpha0 T0 LogML alpha_fit T_fit\r\n'); % Heading of the data file 
    
    % Change to the path containing the behavioural data files
    Dataj = ['./Data/BlackWhiteTest_',num2str(subj+12),'.txt'];
    
    parfor ii=1:size(alpha0_vals,1)
        
        tic; % tic toc to measure run time

        alpha0 = alpha0_vals(ii);
        T0 = T0_vals(ii);

        [alpha_fit,T_fit,LogML] = PWRL_fitting([alpha0,T0],rho,r0,N,Ns,Dataj);
        
        % store the fitted values
        Matj(ii,:) = [alpha0,T0,alpha_fit,T_fit,LogML];
        
        toc
        
    end 
    
    % Record Matj into a text file
    for nn = 1:size(alpha0_vals,1)
        
        % Write trial result to file:   
        fprintf(datafilepointer2,'%1.3f %1.3f %1.4f %1.3f %1.3f %1.3f\r\n', ...
                Matj(nn,1), ...
                Matj(nn,2), ...
                Matj(nn,5), ...
                Matj(nn,3), ...
                Matj(nn,4));

    end
    
    % get the pair (alpha_fit,T_fit) that has the highest logML
    [logML_max, rowsOfMaxes] = max(Matj(:,5));
    
    alpha_bestfit=Matj(rowsOfMaxes, 3); % best fitted value for alpha
    T_bestfit=Matj(rowsOfMaxes, 4); % best fitted value for T
    [aic, bic] = aicbic(logML_max, numOfParams, N); % AIC and BIC scores
    
    MatFittedParams(subj,2)=alpha_bestfit;
    MatFittedParams(subj,3)=T_bestfit;
    MatFittedParams(subj,4)=logML_max;
    MatFittedParams(subj,5) = aic;
    MatFittedParams(subj,6) = bic;
    
    % Write trial result to file:   
    fprintf(datafilepointer1,'%i %1.6f %1.6f %1.2f %1.2f %1.2f\r\n', ...
            subj+12, ...
            alpha_bestfit, ...
            T_bestfit, ...
            logML_max, ...
            aic, ...
            bic); 

    myfilename = ['./Results/FittingFile_subj_',num2str(subj+12),'.mat']; 
    save(myfilename);
    
    % Output message to have a feel of the running progress
    ending_subj = clock;   
    fprintf('****Subject: %d - Elapsed time: %1.2f\r\n',[subj,etime(ending_subj,beg_subj)]);

end

% save fitting results in a '.mat' file
myfilename = './Results/FitTest.mat'; 
save(myfilename,'MatFittedParams');

% Wrapping up
delete(pool)
fclose('all');
ending = clock;
total_duration = etime(ending,beg);



