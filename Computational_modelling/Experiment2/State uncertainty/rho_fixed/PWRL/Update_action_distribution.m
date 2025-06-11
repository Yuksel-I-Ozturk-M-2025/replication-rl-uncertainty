function P=Update_action_distribution(QQ,o,T,rho)
%P_action_state=Update_action_distribution(Q,O(n),T,rho);

% This function update the distribution of taking an action given the
% current state o

P=zeros(2,1);

if o==0
    rho_t=[rho;1-rho]; % probabilities of identifying the current state as dark and light respectively, given that the current state is dark 
elseif o==1
    rho_t=[1-rho;rho]; % probabilities of identifying the current state as dark and light respectively, given that the current state is light 
end

for a=1:2 % (a=2=circle && a=1=square)
    % update the distribution of actions given a state s 
    P(a)= (rho_t(1)/( 1 + exp((QQ(1,3-a)-QQ(1,a))/T) )) + (rho_t(2)/( 1 + exp((QQ(2,3-a)-QQ(2,a))/T) )) ;
end


end


