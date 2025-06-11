function QQ=update_state_action_value(Q,o,a,r,alpha)
%QQ=update_state_action_value(Q,o,a,r,alpha);

% QQ corresponds the state-action value vector at step n+1 (updated)
% Q corresponds the state-action value vectorr at step n (old)

% This function update state-action values

QQ=Q;
% for action a (a=0='dark' response, a=1='light' response; we use a+1 because matrices in matlab cannot start from element 0)
QQ(o+1,a+1) = Q(o+1,a+1) + alpha * (r- Q(o+1,a+1)); % update the value for the current state


end