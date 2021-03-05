function norm_of_a = norma(a,nu)

N = length(a) + 1;

nu_power = nu.^abs(2:N)';
norm_of_a = sum(abs(a).*nu_power);

end