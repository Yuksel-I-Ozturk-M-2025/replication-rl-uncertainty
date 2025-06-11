beg = clock; % starting time

% Add the path to the ifit package folder (MODIFY AS NEEDED)
addpath(genpath('./Computational_modelling/ifit'));

paramVector = zeros(8,1);
paramVector(1) = 200; % Total nb of trial in the experiment - N
paramVector(2) = 1000; % nb of simulations to use to estimate the action likelihoods - Ns
paramVector(3) = 2.5; % reward of correct action - r0
paramVector(4) = 0.65; % level of confidence (in perception) - rho
paramVector(5) = 0.8; % weight placed on the grand average payoff of all past experiences - w
paramVector(6) = 0.8; % probability of choosing the most recent trial (t-1) when we sample the 'b' recent trials - pr
paramVector(7) = 6; % bound of the memory - b
paramVector(8) = 1; % exploration rate - T

[N,Ns,r0,rho,w,pr,b,T] = DispParams(paramVector);
J = 12; % number of subjects
numOfParams = 5; % number of parameters to estimate
numOfParamInit = 30; % Num of repetitions we run the fitting procedure for each model-subject - we randomly generated 30 initial parameter values
NumOfWorkers = 16; % Num of cores to use in the foreach loop

%create parallel pool of workers on the local node
%Ensure that this is the same number as what you requested from the scheduler
pool = parpool('local', NumOfWorkers);

% Define filenames of result file (record results for each parameter value averaged over the number of simulation runs):
datafilename1 = strcat('./Results/FittedParameters.txt'); % name of data file to write to
datafilepointer1 = fopen(datafilename1,'wt'); % open ASCII file for writing
fprintf(datafilepointer1,'Subj rho_fit b_fit w_fit pr_fit T_fit LogML AIC BIC\r\n'); % Heading of the data file 

%%% Initial parameters for the fitting procedure
rng('default');
rng(1); % set seed to be able to replicate the results
w0_vals = rand(numOfParamInit,1); 
pr0_vals = rand(numOfParamInit,1);
T0_vals = rand(numOfParamInit,1);
rho0_vals=(1+rand(numOfParamInit,1))/2; % to make sure rho is initialised between 0.5 and 1 

% Possible values for b
b_vals = (1:10)'; 

%%% Matrix containing fitted_w, fitted_pr and log(ML) values for each starting value of w and pr for one subject and one value of b
N_w0 = length(w0_vals); % number of starting values for alpha
N_pr = length(pr0_vals); % number of starting values for T
N_T = length(T0_vals); % number of starting values for T
N_rho0 = length(rho0_vals); % number of starting values for rho
Matb = zeros(N_w0,9);

% best fitted paramters w and pr for each possible b valu
MatFit_b = zeros(length(b_vals),6);
MatFit_b(:,1) = b_vals; 

% Matrix containing the best fitted values with the corresponding logML for each subject
MatFittedParams = zeros(J, 9);
MatFittedParams(:,1) = (1:J)'; 

beg_subj = clock;

for subj = 1:J 

    % Define filenames of result file (record results for each parameter value averaged over the number of simulation runs):
    datafilename2 = strcat('./Results/FittedParameters_subj',num2str(subj+12),'.txt'); % name of data file to write to
    datafilepointer2 = fopen(datafilename2,'wt'); % open ASCII file for writing
    fprintf(datafilepointer2,'b rho0 w0 pr0 T0 LogML rho_fit w_fit pr_fit T_fit\r\n'); % Heading of the data file  
    
    % Change to the path containing the behavioural data files
    Dataj = ['./Data/BlackWhiteTest_',num2str(subj+12),'.txt'];
    
    for i_b = 1:size(b_vals,1)

        b = b_vals(i_b);

        parfor ii = 1:size(w0_vals,1)
            
            tic; % tic toc to measure run time
            
            w0 = w0_vals(ii);
            pr0 = pr0_vals(ii);            
            T0 = T0_vals(ii); 
            rho0 = rho0_vals(ii);
            [w_fit,pr_fit,T_fit,rho_fit,LogML] = BISAW_fitting([w0,pr0,T0,rho0],b,N,Ns,Dataj);
            
            % store the fitted values
            Matb(ii,:) = [rho0, w0, pr0, T0, rho_fit, w_fit, pr_fit, T_fit, LogML];       

            toc   

        end
        
        % Record Matb into a text file
        for nn = 1:size(w0_vals,1)

            % Write results to file:  
            fprintf(datafilepointer2,'%i %1.3f %1.3f %1.3f %1.3f %1.4f %1.4f %1.4f %1.4f %1.4f\r\n', ...
                    b, ...
                    Matb(nn,1), ...
                    Matb(nn,2), ...
                    Matb(nn,3), ...
                    Matb(nn,4), ...
                    Matb(nn,9), ...
                    Matb(nn,5), ...
                    Matb(nn,6), ...
                    Matb(nn,7), ...
                    Matb(nn,8)); 

        end


        % get the pair (w_fit,pr_fit) that has the highest logML
        [logML_max_b, rowsOfMaxes_b] = max(Matb(:,9));
        
        rho_bestfit_b = Matb(rowsOfMaxes_b, 5); % best fitted value for T
        w_bestfit_b = Matb(rowsOfMaxes_b, 6); % best fitted value for w
        pr_bestfit_b = Matb(rowsOfMaxes_b, 7); % best fitted value for pr
        T_bestfit_b = Matb(rowsOfMaxes_b, 8); % best fitted value for T
        
        MatFit_b(i_b,2) = rho_bestfit_b;
        MatFit_b(i_b,3) = w_bestfit_b;
        MatFit_b(i_b,4) = pr_bestfit_b;
        MatFit_b(i_b,5) = T_bestfit_b;
        MatFit_b(i_b,6) = logML_max_b;

    end

    % get the pair (alpha_fit,T_fit) that has the highest logML
    [logML_max, rowsOfMaxes] = max(MatFit_b(:,6));

    b_bestfit = MatFit_b(rowsOfMaxes,1); % best fitted value for b
    rho_bestfit = MatFit_b(rowsOfMaxes,2); % best fitted value for rho
    w_bestfit = MatFit_b(rowsOfMaxes,3); % best fitted value for w
    pr_bestfit = MatFit_b(rowsOfMaxes,4); % best fitted value for pr
    T_bestfit = MatFit_b(rowsOfMaxes,5); % best fitted value for T
    [aic, bic] = aicbic(logML_max, numOfParams, N); % AIC and BIC scores

    MatFittedParams(subj,2) = rho_bestfit;
    MatFittedParams(subj,3) = b_bestfit;
    MatFittedParams(subj,4) = w_bestfit;
    MatFittedParams(subj,5) = pr_bestfit;
    MatFittedParams(subj,6) = T_bestfit;
    MatFittedParams(subj,7) = logML_max;
    MatFittedParams(subj,8) = aic;
    MatFittedParams(subj,9) = bic;

    % Write trial result to file:   
    fprintf(datafilepointer1,'%i %1.4f %i %1.4f %1.4f %1.4f %1.4f %1.2f %1.2f\r\n', ...
            subj+12, ...
            rho_bestfit, ...
            b_bestfit, ...
            w_bestfit, ...
            pr_bestfit, ...   
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
ending=clock;
total_duration = etime(ending,beg);
