function esv = Compute_one_ESV(Rj, m, w, Prob)
    
    % Numbers of previous payoffs for each action given the current state
    Nj = length(Rj);

    %%%%%%%%%%%%%%%%%%% ESV of a = 0 %%%%%%%%%%%%%%%%%%
    if (Nj == 0)
        
        esv = 0; % no payoff has been recorded then use 0        

    elseif (Nj == 1) 
        
        % Compute ESV
        GrandM0 = mean(Rj);
        esv = (1-w)*GrandM0 + w*GrandM0; 

    else 
        
        Ij = Nj + 1 - sample_discrete(repmat(Prob(:, Nj), [1, m])); % draw a sample of m past trials following the 'Prob' distribution 
                
        % Compute the Grand and sample average payoffs for action a = 0:
        GrandM0 = mean(Rj);
        SampleM0 = mean(Rj(Ij));
        esv = (1-w)*SampleM0 + w*GrandM0; % estimated subjective value of the promoted mapping

    end

end