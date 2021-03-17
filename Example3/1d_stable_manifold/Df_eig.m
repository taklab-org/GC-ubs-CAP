function Df = Df_eig(X)

x = X(1:2); lambda = X(3); v = X(4:5);

%f = [phi_vector_field(x);norm(v,2)^2-1;Dphi_vector_field(x)*v-lambda*v];

Df = [Dphi_vector_field(x) zeros(2,3) ;
      0  0  0  2*v(1)  2*v(2) ;
      D2phi_vector_field(x,v) -v  Dphi_vector_field(x)-lambda*eye(2)];

end
