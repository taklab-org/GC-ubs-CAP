function [rho,C] = exp_decay_a_least_square(a)

abs_a=abs(a);

abs_a = abs_a(abs_a>0);
%abs_a = abs_a(1:m);
m = length(abs_a);

Y=log(abs_a);

k=(1:m)';
sum_k=sum(k);
sum_k2=sum(k.^2);

sum_Y=sum(Y);
sum_kY=sum(k.*Y);
denom=m*sum_k2-sum_k^2;
A=(m*sum_kY-sum_k*sum_Y)/denom;
B=(sum_Y*sum_k2-sum_k*sum_kY)/denom;

C=exp(B);
rho=exp(-A);

end






