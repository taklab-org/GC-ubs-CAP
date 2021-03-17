function Dg = Dg_vector_field(x)

if isintval(x) == 1
    Dg = intval(zeros(3)); Jf = intval(zeros(3));
    % Fix the parameters
    a = intval('0.3'); c = intval('0.7'); delta = intval(9); w = intval('0.02');
else
    Dg = zeros(3); Jf = zeros(3);
    % Fix the parameters
    a = 0.3; c = 0.7; delta = 9; w = 0.02;
end

x0 = x(1); x1 = x(2); x2 = x(3);

f31 = 2*x0*x2 + (1/delta)*(2*c*x0*x2+a*x1*(x0-x1)-x1*(x1-a*x0)+3*w*x0^2);
f32 = (1/delta)*(-(x1-a*x0)*(x0-x1)-x1*(x0-x1)+x1*(x1-a*x0));
f33 = x0^2 + (1/delta)*c*x0^2;

Jf(1,1) = 5*x0^2 - (1-x0^2-x1^2-x2^2);
Jf(1,2) = 2*x0*x1; 
Jf(1,3) = 2*x0*x2; 
Jf(2,1) = 2*x0*(x1+x2);
Jf(2,2) = x0^2; 
Jf(2,3) = x0^2; 
Jf(3,1) = f31;
Jf(3,2) = f32;
Jf(3,3) = f33;

tf0 = x0^3 - (1-x0^2-x1^2-x2^2)*x0;
tf1 = x0^2*x1 + x0^2*x2;
tf2 = x0^2*x2 + (1/delta)*(c*x0^2*x2-x1*(x1-a*x0)*(x0-x1)+w*x0^3);

sum = x0*tf0+x1*tf1+x2*tf2;

Dg(1,1) = Jf(1,1) - sum - x0*(tf0+x0*Jf(1,1)+x1*Jf(2,1)+x2*Jf(3,1));
Dg(1,2) = Jf(1,2) - x0*(tf1+x0*Jf(1,2)+x1*Jf(2,2)+x2*Jf(3,2));
Dg(1,3) = Jf(1,3) - x0*(tf2+x0*Jf(1,3)+x1*Jf(2,3)+x2*Jf(3,3));

Dg(2,1) = Jf(2,1) - x1*(tf0+x0*Jf(1,1)+x1*Jf(2,1)+x2*Jf(3,1));
Dg(2,2) = Jf(2,2) - sum - x1*(tf1+x0*Jf(1,2)+x1*Jf(2,2)+x2*Jf(3,2));
Dg(2,3) = Jf(2,3) - x1*(tf2+x0*Jf(1,3)+x1*Jf(2,3)+x2*Jf(3,3));

Dg(3,1) = Jf(3,1) - x2*(tf0+x0*Jf(1,1)+x1*Jf(2,1)+x2*Jf(3,1));
Dg(3,2) = Jf(3,2) - x2*(tf1+x0*Jf(1,2)+x1*Jf(2,2)+x2*Jf(3,2));
Dg(3,3) = Jf(3,3) - sum - x2*(tf2+x0*Jf(1,3)+x1*Jf(2,3)+x2*Jf(3,3));

end