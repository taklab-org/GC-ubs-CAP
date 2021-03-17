function f = f_eig(X)

x = X(1:2); lambda = X(3); v = X(4:5);

f = [phi_vector_field(x);norm(v,2)-1;Dphi_vector_field(x)*v-lambda*v];

end

