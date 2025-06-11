function [A,R,Acc,Qt] = PWRL(O,N,rho,alpha,T,r0)
%[A,R,Acc,Qt] = PWRL(O,N,rho,alpha,T,r0);

% Simulation of my experiment using the PWRL model (see Larsen et al. (2009) )

% A is a verctor whose elements correspond to the action (or mapping) selected at trial n
%A=1 correspond to the circle (correct action for "light" state) and A=0=square (correct action for "dark" state)       

% R is a vector whose elements correspond to the payoff of the promoted mapping

% Q is is a 2(possible states)*2(actions) matrix whose elements correspond to the state-action values

% Acc_avg is a vector that encodes the average accuracy accross the M simulated subjects in each trial

% 1 <-> 'light';
% 0 <-> 'dark';

A = 99*ones(N,1);
R = zeros(N,1);
Acc = zeros(N,1); % accuracy in each trial
Q = ones(2,2);
Qt = zeros(2,2,200); % Q at each time step

for n = 1:N        

    % the subject can only identify the state as the true state with probability rho (o is the identified state)
    i = sample_discrete([rho 1-rho]', 'prob');
    if (i == 1)
        o = O(n); % identify correcly current observation
    else
        o = 1-O(n); % incorrect identification
    end

    % Update the distribution of actions given a state s
    P_action_state = Update_action_distribution1(Q(o+1,:),T);                   

    % select a motor action based on the calculated distribution
    a = sample_discrete(P_action_state, 'prob')-1;     

    % observe the reward
    r = Reward_function(O(n),a,r0,n,N);

    % update the state-action values       
    Q = update_state_action_value(Q,O(n),a,r,r0,alpha,rho);

    % record reward, the selected action and the state-action values
    R(n) = r;
    A(n) = a;
    Qt(:,:,n) = Q;

end

Acc(R>0) = 1; % accuracy vector in the m-th simulation 

end