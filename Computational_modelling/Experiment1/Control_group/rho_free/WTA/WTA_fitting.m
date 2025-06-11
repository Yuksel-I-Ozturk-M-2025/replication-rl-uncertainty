function [alpha_fit,T_fit,rho_fit,LogML] = WTA_fitting(params_init,N,Ns,myfilename)
%[alpha_fit,T_fit,rho_fit,LogML] = WTA_fitting([0.1,0.1,0.1],N,Ns,myfilename);

data = importdata(myfilename);
MatData = data.data; % extract data 

O = MatData(:,5);
A = MatData(:,6);
R = 100*MatData(:,8);

function NeglogML = MLE_RL(free_params,N,Ns,O,A,R)

    
    alpha = free_params(1);
    T = free_params(2);
    rho = free_params(3);
    
    % sampling identification variables before loop
    rng('default');
    Id = reshape(sample_discrete(repmat([rho;1-rho],[1 N*Ns]), 'prob'),Ns,N);

    An_predicted = zeros(N,Ns); %  sampled action from WTA given previous actions up to n-1 to use to estimate the likelihoods of selected actions    
    P_action_state = zeros(2,1);
                
    for ii = 1:Ns
        
        Q_fit = ones(2,2); % initialisation
        
        for n = 1:N   
            
            % the subject can only identify the state as the true state with probability rho (o is the identified state)
            if (Id(ii,n) == 1)
                o = O(n); % identify correcly current observation
            else
                o = 1-O(n); % incorrect identification
            end
            
            % Update the distribution of actions given a state s            
            P_action_state(1) = 1/( 1 + exp((Q_fit(o+1,2)-Q_fit(o+1,1))/T) );
            P_action_state(2) = 1 - P_action_state(1);

            % select a motor action based on the calculated distribution
            An_predicted(n,ii) = sample_discrete(P_action_state, 'prob')-1;

            % update Q            
            Q_fit(o+1,A(n)+1) = Q_fit(o+1,A(n)+1) + alpha * (R(n)- Q_fit(o+1,A(n)+1)); % update the value for the current state

        end 
    
    end
    
    logML = 0;
    for n = 1:N
        logML = logML + log(mean(An_predicted(n,:)==A(n))); % updating the log likelihood 
    end
    NeglogML = -logML;

end

options = optimset('display','off');
options = optimset(options,'MaxIter',200*numel(params_init));
options = optimset(options,'MaxFunEvals',200*numel(params_init));

f = @(free_params)MLE_RL(free_params,N,Ns,O,A,R);
[params,fval] = fminsimpsa(f, params_init, options, [0.0001 0.0001 0.5001], [1 100 0.9999]);

alpha_fit = params(1);
T_fit = params(2);
rho_fit=params(3);
LogML = -fval;
   
end