function [A,R,Acc,Qt] = SARSA(O,N,rho,alpha,T,r0)
%[A,R,Acc,Qt] = SARSA(O,N,rho,alpha,T,r0);

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

    o = O(n); % Current observation

    % Update the distribution of actions given a state s
    P_action_state = Update_action_distribution(Q(o+1,:),T);                   

    % select an action based on the calculated distribution
    a = sample_discrete(P_action_state, 'prob')-1;     

    % observe the reward
    [r,acc] = Reward_function(O(n),a,r0,rho,n,N);

    % update the state-action values       
    Q = update_state_action_value(Q,o,a,r,alpha);

    % record reward, the selected action and the state-action values
    R(n) = r;
    Acc(n) = acc;
    A(n) = a;
    Qt(:,:,n) = Q;
    
    
end

end