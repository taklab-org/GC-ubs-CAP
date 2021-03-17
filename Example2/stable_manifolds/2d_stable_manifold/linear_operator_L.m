function L = linear_operator_L(lambda,N)

if isintval(lambda(1)) == 1
    L = intval(zeros((N+1)^2));
else
    L = zeros((N+1)^2);
end

n1 = (0:N);
LAMBDA1 = diag(lambda(1)*n1);
Id = eye(N+1);

for n2 = 0:N
    L_n2 = lambda(2)*n2*Id + LAMBDA1;
    L(n1+1+n2*(N+1),n1+1+n2*(N+1)) = L_n2;
end

end

