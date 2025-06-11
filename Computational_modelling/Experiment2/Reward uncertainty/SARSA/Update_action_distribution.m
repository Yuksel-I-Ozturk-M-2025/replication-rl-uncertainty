function P = Update_action_distribution(QQ,T)
%P_action_state = Update_action_distribution(QQ,T);

% This function update the distribution of motor and gating action given a state s

P = zeros(2,1);

for a = 1:2 % (a=2=circle && a=1=square)
    % update the distribution of actions given a state s 
    P(a) = 1/( 1 + exp((QQ(3-a)-QQ(a))/T) );
end

end


