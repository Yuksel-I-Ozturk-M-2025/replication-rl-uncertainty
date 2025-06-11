function [r, acc] = Reward_function(o, a, r0, n, N)
%r=Reward_function(O(n),a,r0,n,N);

% this function gives the observed reward at step n+1

% r is the observed reward
% o is the observation at step n
% a is the selected motor action 

acc = 0;

if (n<N/2+1) % Before the switch

    % observe the reward
    if (o == a)
        r = r0; 
        acc = 1;
    else
        r = 0; 
    end   
    
else % After the switch

    % observe the reward
    if (o == a)
        r = 0; 
    else
        r = r0;
        acc = 1;
    end  
    
end

end