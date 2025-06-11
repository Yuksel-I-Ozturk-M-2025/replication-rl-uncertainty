function [alpha_fit,T_fit,LogML] = PWRL_fitting(params_init,rho,r0,N,Ns,myfilename)
%[alpha_fit,T_fit,LogML] = PWRL_fitting([0.1,0.1],rho,N,Ns,'Simulation1.mat');

data = importdata(myfilename);
MatData = data.data; % extract data 

O = MatData(:,5);
A = MatData(:,6);
A(101:200) = 1-A(101:200); % Because A (i.e. resp) encodes whether it is the dark (A=0) or light response (A=1), and this encoding remains the same after the switch (so needed to adjust A) 
R = MatData(:,8);

function NeglogML = MLE_RL(free_params,N,Ns,rho,r0,O,A,R)
    
    % set the seed
    rng('default');
    
    alpha = free_params(1);
    T = free_params(2);    
    
    % sampled action from WTA given previous actions up to n-1 to use to estimate the likelihoods of selected actions  
    An_predicted = zeros(N,Ns); 
    
    for ii = 1:Ns
        
        Q_fit = ones(2,2); % initialisation
        
        for n = 1:N   
            
            % Update the distribution of actions given a state s            
            P_action_state = Update_action_distribution(Q_fit,O(n),T,rho);

            % select a motor action based on the calculated distribution
            An_predicted(n,ii) = sample_discrete(P_action_state, 'prob')-1;

            % update Q            
            Q_fit = update_state_action_value(Q_fit,O(n),A(n),R(n),r0,alpha,rho);

        end 
    
    end
    
    % Compute -loglik 
    B = (An_predicted == A);
    NeglogML = -sum(log(mean(B, 2)));

end

options = optimset('display','off');
options = optimset(options,'MaxIter',200*numel(params_init));
options = optimset(options,'MaxFunEvals',200*numel(params_init));

f = @(free_params)MLE_RL(free_params,N,Ns,rho,r0,O,A,R);  
[params,fval] = fminsimpsa(f, params_init, options, [0.0001 0.0001],[1 100]);

alpha_fit=params(1);
T_fit=params(2);
LogML=-fval;

end