function [r, acc] = Reward_function(o, a, r0)
%[r, acc]=Reward_function(O(n),a,r0);

% this function gives the observed reward at step n+1

% r is the observed reward
% o is the observation at step n
% a is the selected motor action 

acc = 0;

% observe the reward
if (o == a)
    r = r0; 
    acc = 1;
else
    r = 0; 
end  

end