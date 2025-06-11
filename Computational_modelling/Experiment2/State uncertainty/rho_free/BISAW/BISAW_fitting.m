function [w_fit, pr_fit, T_fit, rho_fit, LogML] = BISAW_fitting(params_init,b,N,Ns,myfilename)
%[w_fit,pr_fit,T_fit,rho_fit,LogML] = BISAW_fitting([0.8,0.8,1,0.7],b,N,Ns,myfilename);

data = importdata(myfilename);
MatData = data.data; % extract data 

O = MatData(:,5);
A = MatData(:,6);
A(101:200) = 1-A(101:200);
R = MatData(:,8); 

% Vector containing trial indices
Trials = (1:N)';

function NeglogML = MLE_RL(free_params,b,N,Ns,O,A,R)
    
    rng('default');
   
    w = free_params(1); % weight placed on the grand average payoff of all past experiences
    pr = free_params(2); % probability of choosing the most recent trial (t-1) when we sample the 'mu' recent trials
    T = free_params(3); % exploration rate
    rho = free_params(4); % uncertainty level 
    
    ESV_fit = zeros(2,1);
    An_predicted = zeros(N,Ns); %  sampled action from BI-SAW given previous actions up to n-1 to use to estimate the likelihoods of selected actions
    
    % constructing the distribution that will serve in drawing an m-sample with the most recent trial (n-1) having a proba pr to
    % be selected whereas each other trial has a proba (1-pr)/n-2 of being chosen
    % We initialise it outside the loop to speed up the program
    P0 = (ones(b,1)-pr)/(b-1); 
    P0(1) = pr;
    Prob = repmat(P0,1,N);
    if b>1
        for n = 2:b-1
            Prob(:, n) = 0;
            Prob(2:n, n) = (ones(n-1,1)-pr)/(n-1); 
            Prob(1,n) = pr;
        end
    end
    Prob(:,1) = 0;
    Prob(1,1) = 1;

    % Initialise probability distribution of actions given a known state
    P_action_state = zeros(2,1); 
    
    % sampling identification variables before loop (note that I sample for all Ns responses so there will be still ramdomness)
    Id = reshape(sample_discrete(repmat([rho;1-rho],[1 N*Ns])), Ns, N);
    
    % Initialise vector of identified states 
    Id_states = zeros(N,1); 
                
    for ii = 1:Ns
        
        for n = 1:N   
            
            % the subject can only identify the state as the true state with probability rho (o is the identified state)
            if (Id(ii, n) == 1)
                o = O(n); % identify correcly current observation
            else
                o = 1-O(n); % incorrect identification
            end 
            Id_states(n) = o;
            
            % Reward history for each action
            R0 = R((Trials < n) & (A==0) & (Id_states==o)); % reward history for action a = 0 and current state O(n)
            R1 = R((Trials < n) & (A==1) & (Id_states==o)); % reward history for action a = 1 and current state O(n)
            
            % Numbers of previous payoffs for each action given the current state
            N0 = length(R0);
            N1 = length(R1);

            %%%%%%%%%%%%%%%%%%% ESV of a = 0 %%%%%%%%%%%%%%%%%%
            if (N0 == 0)

                ESV_fit(1) = 0; % no payoff has been recorded then use 0        

            elseif (N0 == 1) 

                % Compute ESV
                ESV_fit(1) = R0(1); 

            else 

                % draw a sample of m past trials following the 'Prob' distribution 
                I0 = N0 + 1 - sample_discrete(repmat(Prob(:, N0), [1, b])); 

                % Compute the Grand and sample average payoffs for action a = 0:
                GrandM0 = mean(R0);
                SampleM0 = mean(R0(I0));
                ESV_fit(1) = (1-w)*SampleM0 + w*GrandM0; % estimated subjective value of the promoted mapping

            end
            
            %%%%%%%%%%%%%%%%%%% ESV of a = 1 %%%%%%%%%%%%%%%%%%
            if (N1 == 0)

                ESV_fit(2) = 0; % no payoff has been recorded then use 0        

            elseif (N1 == 1) 

                % Compute ESV
                ESV_fit(2) = R1(1); 

            else 

                % draw a sample of m past trials following the 'Prob' distribution 
                I1 = N1 + 1 - sample_discrete(repmat(Prob(:, N1), [1, b]));

                % Compute the Grand and sample average payoffs for action a = 1:
                GrandM1 = mean(R1);
                SampleM1 = mean(R1(I1));
                ESV_fit(2) = (1-w)*SampleM1 + w*GrandM1; % estimated subjective value of the promoted mapping

            end
            

            % Update the distribution of actions 
            for a = 1:2 % (a=2=circle && a=1=square)
                P_action_state(a) = 1/( 1 + exp((ESV_fit(3-a)-ESV_fit(a))/T) );
            end

            % select an action based on the calculated distribution
            An_predicted(n, ii) = sample_discrete(P_action_state)-1; 
             
        end 
    
    end
    
    % Compute -loglik 
    B = (An_predicted == A);
    NeglogML = -sum(log(mean(B, 2)));
    
end

options = optimset('display','off');
options = optimset(options,'MaxIter',200*numel(params_init));
options = optimset(options,'MaxFunEvals',200*numel(params_init));

f = @(free_params)MLE_RL(free_params,b,N,Ns,O,A,R);
[params,fval] = fminsimpsa(f, params_init, options, [0.0001 0.0001 0.0001 0.5001], [1 1 100 0.9999]);

w_fit = params(1); % weight placed on the grand average payoff of all past experiences
pr_fit = params(2); % probability of choosing the most recent trial (t-1) when we sample the 'mu' recent trials
T_fit = params(3); % exploration rate
rho_fit = params(4); % uncertainty level

LogML = -fval;

end