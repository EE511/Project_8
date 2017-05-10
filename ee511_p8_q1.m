%%--------------------------------------------------------------------------
%%Project-8:: Question - 1
%%To simulate Monte Carlo stratification and importance sampling
%%Author                Date               Revision
%%Rajasekar Raja     05/07/17         Initial Revision
%%--------------------------------------------------------------------------
function [ ] = monte_carlo_simulation(No_of_samples)
% Estimating mean and variance using simple Monte Carlo
g = @(x1,x2)exp(5.*abs(x1-5) + 5.*abs(x2-5));
a = -1;
b = 1;
r_number1 = rand(1,No_of_samples);
r_number2 = rand(1,No_of_samples);
X = g(r_number1,r_number2);
display(['Mean ',num2str(mean(X)),' and Variance ',num2str(2*std(X)/sqrt(No_of_samples)),' using simple Monte Carlo']);

%Stratified sampling
K = 20; 
Nij = No_of_samples/K;
for i = 1:K
    for j = 1:K
        XS = g((i-1+rand(1,Nij))/K,(j-1+rand(1,Nij))/K);
        XSb(i,j) = mean(XS); 
        SS(i,j) = var(XS);
    end
end
SST = mean(mean(SS/No_of_samples));
display(['Mean ',num2str(mean(mean(XSb))),' and Variance ',num2str(2*sqrt(SST)),' using Stratified Sampling']);

% Estimating mean and variance using Importance Sampling
e = exp(1);% Importance sampling
u3 = rand(1,No_of_samples);
u4 = rand(1,No_of_samples);
X1 = log(1+(e-1)*u3);
X2 = log(1+(e-1)*u4);
T = (e-1)^2*exp(5.*abs(X1-5) + 5.*abs(X2-5)-(X1+X2));
display(['Mean ',num2str(mean(T)),' and Variance ',num2str(2*std(T)/sqrt(No_of_samples)),' using Importance Sampling']);
pl = 0;
pu = 1;
ql = 0;
qu = 1;
fun = @(p,q) exp(5.*abs(p-5) + 5.*abs(q-5));                  % theoretical integration of funtion(x,y)
theo_val = integral2 (@(p,q)fun(p,q),pl,pu,ql,qu);
disp(['Theoretical integral value ',num2str(theo_val)]);

%%%%%%%%%%%%%%% Part b %%%%%%%%%%%%%%%

No_of_samples = 1000; 
a = -1; 
b = 1; 
X12 = rand(1,No_of_samples); X22=rand(1,No_of_samples); 
V=(b-a)*(b-a)*cos(pi + 5*(a + (b-a)*X12) + 5*(a + (b-a)*X22));
display(['Mean ',num2str(mean(V)),' and Variance ',num2str(2*std(V)/sqrt(No_of_samples)),' using simple Monte Carlo']);
g = @(x1,x2)cos(pi + 5.*x1 + 5.*x2);

K = 20; Nij = No_of_samples/K;                       % With Stratified sampling
for i = 1:K
    for j = 1:K
        r1 = (b-a).*rand(1,Nij) + a;
        r2 = (b-a).*rand(1,Nij) + a;
        XS = g((i-1+r1)/K,(j-1+r2)/K);
        XSb(i,j) = mean(XS); 
        SS(i,j) = var(XS);
    end
end
SST = mean(mean(SS/No_of_samples));
display(['Mean ',num2str(mean(mean(XSb))),' and Variance ',num2str(2*sqrt(SST)),' using simple Monte Carlo']);

% Estimating mean and variance using Importance Sampling

e = exp(1);                               % Importance sampling
r3 = (b-a).*rand(1,Nij) + a;
r4 = (b-a).*rand(1,Nij) + a;
X12 = log(1+(e-1)*r3);
X22 = log(1+(e-1)*r4);
T = (e-1)^2*cos(pi + 5.*X12 + 5.*X22-(X12+X22));
display(['Mean ',num2str(mean(T)),' and Variance ',num2str(2*std(T)/sqrt(No_of_samples)),' using simple Monte Carlo']);

pl = -1;
pu = 1;
ql = -1;
qu = 1;
fun = @(p,q) cos(pi + 5*p + 5*q);                  % theoretical integration of funtion(x,y)
theo_val = integral2 (@(p,q)fun(p,q),pl,pu,ql,qu);
disp(['Theoretical integral value ',num2str(theo_val)]);
