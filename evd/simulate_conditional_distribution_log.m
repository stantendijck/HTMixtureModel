function y = simulate_conditional_distribution_log(n,alpha)

y = -alpha*log(rand(n,1).^(1/(alpha-1))-1);

end