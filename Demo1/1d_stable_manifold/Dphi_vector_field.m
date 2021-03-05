function Dphi = Dphi_vector_field(x,par)

Dphi = zeros(2);

x1 = x(1); x2 = x(2);

rho1 = par(2); 
rho2 = par(3); 
c = par(4);
c1 = par(5);
c2 = par(6);

Dphi(1,1) = 3*x1^2 - 2*(rho1+rho2)*x1 + rho1*rho2 - 3*c*x1^2*x2 - 2*c1*x1*x2;
Dphi(1,2) = - c*x1^3 - c1*x1^2;
Dphi(2,1) = -x1*x2 + 2*c*x1*x2^2 + 2*c2*x1*x2^3;
Dphi(2,2) = -(1/2)*x1^2 + (1/2)*rho1*rho2 + 2*c*x1^2*x2 + 3*c2*x1^2*x2^2;

end

