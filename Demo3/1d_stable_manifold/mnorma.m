function matrix_norm=mnorma(M,nu)

N = length(M) + 1 ;

nu_power = nu.^abs(2:N);

matrix_norm = intval(max(mag(sum(abs(M).*repmat(nu_power,N-1,1)')./nu_power)));

end