function r=Reward_function(o,a,r0,n,N)
%r=Reward_function(O(n),a,r0,n,N);

% this function gives the observed reward at step n+1

% r is the observed reward
% o is the observation at step n
% a is the selected motor action 

if (n<N/2+1) % Before the switch

    % observe the reward
    if (o==1 && a==1)
        r=r0; % positive reward if "light" response (for e.g., circle) chosen with "light" state
    elseif (o==1 && a==0)
        r=0; % no reward if "dark" response (for e.g., square) chosen with "light" state 
    elseif (o==0 && a==0)
        r=r0;
    else
        r=0;
    end   
    
else % After the switch

    % observe the reward
    if (o==1 && a==0)
        r=r0; % positive reward if "dark" response (for e.g., circle) chosen with "light" state
    elseif (o==1 && a==1)
        r=0; % no reward if "light" response (for e.g., square) chosen with "light" state 
    elseif (o==0 && a==1)
        r=r0;
    else
        r=0;
    end  
    
end

end