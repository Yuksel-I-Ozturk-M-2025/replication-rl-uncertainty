beg = clock;

tic; % starting time

paramVector = zeros(8,1);
paramVector(1) = 200; % N - total nb of steps in the simulation
paramVector(2) = 100; % Ns - nb of simulations to use to estimate the action likelihoods - Ns = 100
paramVector(3) = 10; % r0 - reward of correct action - r0 = 10
paramVector(4) = 0.7; % rho - level of confidence (in perception)
paramVector(5) = 0.8; % w - weight placed on the grand average payoff of all past experiences - w = 0.8
paramVector(6) = 0.8; % pr - probability of choosing the most recent trial (t-1) when we sample the 'mu' recent trials - pr = 0.8
paramVector(7) = 6; % b - bound of the memory - b = 6
paramVector(8) = 1; % exploration rate - T

[N, Ns, r0, rho, w, pr, b, T] = DispParams(paramVector);
J = 31; % number of subjects
nb_rep = 100; % Number of simulated learning experiences

% Define filenames of result file (record results for each parameter value averaged over the number of simulation runs):
datafilename1  =  strcat('./Results/Avg_performance.txt'); % name of data file to write to
datafilepointer1  =  fopen(datafilename1,'wt'); % open ASCII file for writing
fprintf(datafilepointer1,'N rho b w pr T avg_accuracy relative_reward\r\n'); % Heading of the data file 

% initialisation
reward_rate = 0; 
% Proportion of accurate responses
Acc_avg = zeros(nb_rep, 1); % encodes the overall accuracy for each participant
Acc_BISAW = zeros(N, 31); % For each subject's fitted parameters, one set of learning curves was obtained by averaging over nb_rep (100) simulations of the model 
Acc_all = zeros(N, nb_rep); % Acc in nb_rep (100) simulations for each subject
% Proportion of rewarded responses
R_rel_avg = zeros(nb_rep, 1); % encodes the overall relative reward for each participant 
R_rel_BISAW = zeros(N, 31); % For each subject's fitted parameters, one set of learning curves was obtained by averaging over nb_rep (100) simulations of the model 
R_rel_all = zeros(N, nb_rep); % relative reward in nb_rep (100) simulations for each subject


%%% set of parameters to use
load('./Results/FitTest.mat', 'MatFittedParams');
b_vals = MatFittedParams(:,2);
w_vals = MatFittedParams(:,3);
pr_vals = MatFittedParams(:,4);
T_vals = MatFittedParams(:,5);

for subj = 1:J

    b = b_vals(subj);
    w = w_vals(subj);
    pr = pr_vals(subj);
    T = T_vals(subj);
    
    % extract observations for the current subj
    Dataj = ['./Data/StateUncertaintyPhase2_',num2str(subj+2),'.txt'];
    data = importdata(Dataj);
    MatData = data.data; % extract data 
    O = MatData(:,5);

    for ii = 1:nb_rep

        [A, R, Acc, ESV] = BISAW(O, N, b, pr, b, w, rho, r0, T);
        % Overall accuracy
        Acc_avg(ii) = mean(Acc); 
        Acc_all(:, ii) = Acc;
        % Overall relative reward
        R_rel_avg(ii) = mean(R/r0); 
        R_rel_all(:, ii) = R/r0;

    end

    Acc_BISAW(:, subj) = mean(Acc_all, 2);
    R_rel_BISAW(:, subj) = mean(R_rel_all, 2);

    % Write trial result to file:   
    fprintf(datafilepointer1,'%i %1.2f %i %1.4f %1.4f %1.4f %1.4f %1.4f\r\n', ...
            N, ...
            rho, ...     
            b, ...
            w, ...
            pr, ...
            T, ...
            mean(Acc_avg), ...
            mean(R_rel_avg)); 

    toc   

end

% save fitting results in a '.mat' file
save('./Results/Perf_BISAW_fit.mat','Acc_BISAW','R_rel_BISAW');
    

% Wrap up
fclose('all');
ending = clock;
total_duration = etime(ending,beg);



