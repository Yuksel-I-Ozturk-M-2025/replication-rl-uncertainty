beg = clock;

tic; % starting time

paramVector = zeros(5,1);
paramVector(1) = 200; % N - total nb of steps in the simulation
paramVector(2) = 0.65; % rho - level of confidence (in perception)
paramVector(3) = 0.02; % alpha - learning rate
paramVector(4) = 0.05; % T - temperature
paramVector(5) = 2.5; % r0 - reward value for correct responses

[N,rho,alpha,T,r0] = DispParams(paramVector);
J = 12; % number of subjects
nb_rep = 100; % Number of simulated learning experiences

% Define filenames of result file (record results for each parameter value averaged over the number of simulation runs):
datafilename1  =  strcat('.\Results\Avg_performance.txt'); % name of data file to write to
datafilepointer1  =  fopen(datafilename1,'wt'); % open ASCII file for writing
fprintf(datafilepointer1,'N rho alpha T avg_accuracy relative_reward\n'); % Heading of the data file 

% initialisation
reward_rate = 0; 
% Proportion of accurate responses
Acc_avg = zeros(nb_rep, 1); % encodes the overall accuracy for each participant
Acc_PWRL = zeros(N, J); % For each subject's fitted parameters, one set of learning curves was obtained by averaging over nb_rep (100) simulations of the model 
Acc_all = zeros(N, nb_rep); % Acc in nb_rep (100) simulations for each subject
% Proportion of rewarded responses
R_rel_avg = zeros(nb_rep, 1); % encodes the overall relative reward for each participant 
R_rel_PWRL = zeros(N, J); % For each subject's fitted parameters, one set of learning curves was obtained by averaging over nb_rep (100) simulations of the model 
R_rel_all = zeros(N, nb_rep); % relative reward in nb_rep (100) simulations for each subject

%%% set of parameters to use
load('.\Results\FitControl.mat', 'MatFittedParams');
alpha_vals = MatFittedParams(:, 2);
T_vals = MatFittedParams(:, 3);

for subj = 1:J

    alpha = alpha_vals(subj);
    T = T_vals(subj);
    
    % extract observations for the current subj
    Dataj = ['.\Data\BlackWhiteTest_',num2str(subj),'.txt'];
    data = importdata(Dataj);
    MatData = data.data; % extract data 
    O = MatData(:,5);

    for ii = 1:nb_rep

        [A, R, Acc, Qt] = PWRL(O, N, rho, alpha, T, r0);
        % Overall accuracy
        Acc_avg(ii) = mean(Acc); 
        Acc_all(:, ii) = Acc;
        % Overall relative reward
        R_rel_avg(ii) = mean(R/r0); 
        R_rel_all(:, ii) = R/r0;

    end

    Acc_PWRL(:, subj) = mean(Acc_all, 2);
    R_rel_PWRL(:, subj) = mean(R_rel_all, 2);

    % Write trial result to file:   
    fprintf(datafilepointer1,'%i %1.2f %1.4f %1.4f %1.4f %1.4f\n', ...
            N, ...
            rho, ...
            alpha, ...
            T, ...
            mean(Acc_avg), ...
            mean(R_rel_avg)); 

    toc   

end

% save fitting results in a '.mat' file
save('.\Results\Perf_PWRL_fit.mat','Acc_PWRL','R_rel_PWRL');
    


fclose('all');

ending = clock;

total_duration = etime(ending,beg);



