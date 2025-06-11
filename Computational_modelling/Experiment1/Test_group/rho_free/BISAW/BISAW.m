function [A,R,Acc,ESVt] = BISAW(O,N,m,pr,b,w,rho,r0,T)
%[A,R,Acc,ESVt] = BISAW(O,N,m,pr,b,w,rho,r0,T);

% Simulation of my experiment using BI-SAW model (see Chen et al. (2011))

% A is a verctor whose elements correspond to the selected action on trial n
%A=1 corresponds to the correct action for "light" state (O = 1) before the switch
%A=0 corresponds to the correct action for "dark" state (O = 0) before the switch  

% For Observations (O):
%1 <-> 'light';
%0 <-> 'dark';

% R is the vectore whose elements correspond to the received payoffs

% ESV is a matrix whose elements correspond to the estimated subjective values of the state-action pairs on each trial. 
% ESV(s,a,n): estimated subjective value for the state-action pair (s,a) at time step n

% Ij is a vector whose elements correspond to the drawn trials that serve to calculate the average payoff for each state-action pair.
   
A = 99*ones(N,1);
R = zeros(N,1);
ESV = zeros(2,1);
ESVt = zeros(2, N);
Acc = zeros(N,1); % accuracy in each trial
P_action_state = zeros(2,1); % Probability distribution of actions given a known state
Id_states = zeros(N,1); % vector of identified states 

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

% Vector containing trial indices
Trials = (1:N)';

for n = 1:N

    % the subject can only identify the state as the true state with probability rho (o is the identified state)
    i = sample_discrete([rho; 1-rho]);
    if (i==1)
        o = O(n); % identify correcly current observation
    else
        o = 1-O(n); % incorrect identification
    end 
    Id_states(n) = o;
    
    % Reward history for each action
    R0 = R((Trials < n) & (A==0) & (Id_states==o)); % reward history for action a = 0 and current state O(n)
    R1 = R((Trials < n) & (A==1) & (Id_states==o)); % reward history for action a = 1 and current state O(n)

    %%% Compute the ESV of each action
    ESV(1) = Compute_one_ESV(R0, m, w, Prob); % estimated subjective value of the promoted mapping
    ESV(2) = Compute_one_ESV(R1, m, w, Prob); % estimated subjective value of the promoted mapping
    
    % Update the distribution of actions 
    for a = 1:2 % (a=2=circle && a=1=square)
        P_action_state(a) = 1/( 1 + exp((ESV(3-a)-ESV(a))/T) );
    end
    
    % select an action based on the calculated distribution
    a = sample_discrete(P_action_state)-1; 
    
    % observe the reward
    [r, acc] = Reward_function(O(n), a, r0, n, N);
    
    % record reward, the selected action and the state-action values
    R(n) = r;
    Acc(n) = acc;
    A(n) = a;
    ESVt(:,n) = ESV;

end

end