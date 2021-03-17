function phi = phi_vector_field(x,par)

phi = zeros(2,1);

x1 = x(1); x2 = x(2);

rho1 = par(2); 
rho2 = par(3); 
c = par(4);
c1 = par(5);
c2 = par(6);

phi(1) = x1^3 - (rho1+rho2)*x1^2 + rho1*rho2*x1 - c*x1^3*x2 - c1*x1^2*x2;
phi(2) = -(1/2)*x1^2*x2 + (1/2)*rho1*rho2*x2 + c*x1^2*x2^2 + c2*x1^2*x2^3;

end