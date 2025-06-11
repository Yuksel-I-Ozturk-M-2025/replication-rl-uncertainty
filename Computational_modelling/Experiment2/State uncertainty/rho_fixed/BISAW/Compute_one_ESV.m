function esv = Compute_one_ESV(Rj, b, m, w, pr, prob_lastb)
    
    % Numbers of previous payoffs for each action given the current state
    Nj = length(Rj);

    %%%%%%%%%%%%%%%%%%% ESV of a = 0 %%%%%%%%%%%%%%%%%%
    if (Nj == 0)
        
        esv = 0; % no payoff has been recorded then use 0        

    elseif (Nj == 1) 
              
        % Compute the Grand and sample average payoffs
        GrandM0 = mean(Rj);
        
        % Compute ESV
        esv = (1-w)*GrandM0 + w*GrandM0; 

    else 
        
        % Prepare the probabilities to construct the small sample of past experience
        Prob0 = zeros(Nj, 1);
        if (Nj < b) 
            % constructing the distribution that will serve in drawing an m-sample with the most recent trial (n-1) having a proba rho to
            % be selected whereas the other trial having a proba (1-rho)/n-2 of being chosen
            Prob0 = (ones(Nj, 1) - pr)/(Nj-1); 
            Prob0(Nj) = pr;
        else
            % constructing the distribution that will serve in drawing an m-sample with the most recent trial (n-1) having a proba rho to
            % be selected whereas the other trial having a proba (1-rho)/b-1 of being chosen
            Prob0((Nj-b+1):Nj) = prob_lastb;
            Prob0(Nj) = pr;
        end

        I0 = sample_discrete(repmat(Prob0, [1 m])); % draw a sample of m past trials following the 'Prob' distribution
                
        % Compute the Grand and sample average payoffs for action a = 0:
        GrandM0 = mean(Rj);
        SampleM0 = mean(Rj(I0));
        esv = (1-w)*SampleM0 + w*GrandM0; % estimated subjective value of the promoted mapping

    end

end