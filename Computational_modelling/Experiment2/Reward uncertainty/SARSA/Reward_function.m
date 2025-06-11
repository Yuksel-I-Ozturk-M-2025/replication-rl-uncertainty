function [r,acc] = Reward_function(o,a,r0,rho,n,N)
%[r,acc] = Reward_function(O(n),a,r0,rho,n,N);

% this function gives the observed reward at step n+1

% r is the observed reward
% acc is the accuracy in step n
% o is the observation in step n
% a is the selected motor action 

acc = 0;
if (n<N/2+1) % Before the switch

    % observe the reward
    if (o==1 && a==1)
        r = r0*binornd(1,rho); % positive reward with freq rho (65%) if "light" response (for e.g., circle) chosen with "light" state
        acc = 1;
    elseif (o==1 && a==0)
        r = r0*binornd(1,1-rho); % no reward with freq rho (65%) if "dark" response (for e.g., square) chosen with "light" state 
    elseif (o==0 && a==0)
        r = r0*binornd(1,rho);
        acc = 1;
    else
        r = r0*binornd(1,1-rho);
    end   
    
else % After the switch

    % observe the reward
    if (o==1 && a==0)
        r = r0*binornd(1,rho); % positive reward with freq rho (65%) if "dark" response (for e.g., circle) chosen with "light" state
        acc = 1;
    elseif (o==1 && a==1)
        r = r0*binornd(1,1-rho); % no reward with freq rho (65%) if "light" response (for e.g., square) chosen with "light" state 
    elseif (o==0 && a==1)
        r = r0*binornd(1,rho);
        acc = 1;
    else
        r = r0*binornd(1,1-rho);
    end  
    
end

end