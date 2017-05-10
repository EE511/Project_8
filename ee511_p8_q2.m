%%--------------------------------------------------------------------------
%%Project-8:: Question - 2
%%To simulate Gibbs Sampling
%%Author                Date               Revision
%%Rajasekar Raja     05/07/17         Initial Revision
%%--------------------------------------------------------------------------

function [ ] = gibbs_sampling(No_of_samples)
num_rand_var = 3; 
rand_var_vec = 4*ones(1,num_rand_var);
given_mean=1;
for part = 1:2
    if (part ==1)
        constant = 15;
    else
        constant = 1;
    end
    for iter = 1 : No_of_samples 
        rand_var_index = ceil(num_rand_var*rand);   
        S = sum(rand_var_vec) - rand_var_index*rand_var_vec(rand_var_index);    
        if (part ==1)
            rand_var_vec(rand_var_index) = max(constant-S,0)-log(rand)/given_mean;    
        else
            rand_var_vec(rand_var_index) = log(rand)/given_mean- max(constant-S,0);
        end
        sum_rand_var(iter) = S + rand_var_vec(rand_var_index);
    end
    if (part ==1)
        disp(['E[X1+2X2+3X3|X1+2X2+3X3>15]: The Mean of the estimate is ',num2str(mean(sum_rand_var))]);
    else
        disp(['E[X1+2X2+3X3|X1+2X2+3X3<1]: The Mean of the estimate is ',num2str(mean(sum_rand_var))]);
    end
end
end