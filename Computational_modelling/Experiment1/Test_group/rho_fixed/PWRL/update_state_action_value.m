function QQ=update_state_action_value(Q,o,a,r,r0,alpha,rho)
%QQ=update_state_action_value(Q,O(n),a,r,r0,alpha,rho);

% QQ corresponds the state-action value vector at step n+1 (updated)
% Q corresponds the state-action value vectorr at step n (old)

% This function update state-action values

QQ=Q;

if o==0
    rho_t=[rho;1-rho]; % probabilities of identifying the current state as dark and light respectively, given that the current state is dark 
elseif o==1
    rho_t=[1-rho;rho]; % probabilities of identifying the current state as dark and light respectively, given that the current state is light 
end

if (r==0)
    
    for oo=1:2
        % for action a (a=0='dark' response, a=1='light' response; we use a+1 because matrices in matlab cannot start from element 0)
        QQ(oo,a+1) = Q(oo,a+1) + alpha * ( ( (r0- Q(oo,a+1))*rho_t(oo) ) / ( sum( (r0-Q(:,a+1)).*rho_t ) ) ) * (0- Q(oo,a+1)); % update the value for the current state
    end 

elseif (r==r0)
    
    for oo=1:2
        % for action a (a=0='dark' response, a=1='light' response; we use a+1 because matrices in matlab cannot start from element 0)
        QQ(oo,a+1) = Q(oo,a+1) + alpha * ( ( Q(oo,a+1)*rho_t(oo) ) / ( sum( Q(:,a+1).*rho_t ) ) ) * (r0- Q(oo,a+1)); % update the value for the current state
    end
       
end


end