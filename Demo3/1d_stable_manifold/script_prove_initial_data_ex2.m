
clear

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Proof of the fixed point %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

bx = [0.8861081289780320;0.6192579489210105];

ibx = intval(bx);

phi = phi_vector_field(ibx);

Dphi = Dphi_vector_field(ibx);

A = inv(Dphi);

Y0 = norm(A*phi,inf);

disp(['Y0 = ',num2str(sup(Y0))])

r_star = 1e-6;

h = infsup(-r_star,r_star);

b1 = bx(1) + h;
b2 = bx(2) + h;

f1_x1x1 = 1/2 - (15*b1^4)/2 + b2 + b1^2*b2 - 30*b1^4*b2 + 45*b1^8*b2 + (3*b2^2)/2 - 2*b2^3 + 30*b1^4*b2^3 + b2^5;
f1_x1x2 = b1 + b1^3/3 - 6*b1^5 + 5*b1^9 + 3*b1*b2 - 6*b1*b2^2 + 18*b1^5*b2^2 + 5*b1*b2^4;
f1_x2x2 = (3*b1^2)/2 - (9*b2)/2 - 6*b1^2*b2 + 6*b1^6*b2 + 10*b1^2*b2^3;
f2_x1x1 = b1/2 - 5*b1^3 + (21*b1^5)/2 + 90*b1^7 - 117*b1^11 - 40*b1^3*b2 + (23*b1*b2^2)/2 + 10*b1^3*b2^2 - 90*b1^7*b2^2 - 5*b1^3*b2^4;
f2_x1x2 = -10*b1^4 + (3*b2)/2 + (23*b1^2*b2)/2 + 5*b1^4*b2 - (45*b1^8*b2)/2 - 3*b2^3 - 5*b1^4*b2^3 + (3*b2^5)/2;
f2_x2x2 = (3*b1)/2 + (23*b1^3)/6 + b1^5 - (5*b1^9)/2 - 9*b1*b2^2 - 3*b1^5*b2^2 + (15*b1*b2^4)/2;

A11 = A(1,1); A12 = A(1,2); A21 = A(2,1); A22 = A(2,2); 

Z2_1 = norm(A11*f1_x1x1 + A12*f2_x1x1,inf) ...
     + 2*norm(A11*f1_x1x2 + A12*f2_x1x2,inf) ...
       + norm(A11*f1_x2x2 + A12*f2_x2x2,inf);

Z2_2 =   norm(A21*f1_x1x1 + A22*f2_x1x1,inf) ...
     + 2*norm(A21*f1_x1x2 + A22*f2_x1x2,inf) ...
       + norm(A21*f1_x2x2 + A22*f2_x2x2,inf);

Z2 = max([sup(Z2_1) sup(Z2_2)]);

disp(['Z2 = ',num2str(Z2)])

Y0 = intval(Y0); Z0 = intval(0); Z1 = intval(0); Z2 = intval(Z2);

if inf(1-Z0-Z1)>0 
  if inf((1-Z0-Z1)^2-4*Y0*Z2) > 0  
    rmin=sup(((1-Z0-Z1) - sqrt((1-Z0-Z1)^2-4*Y0*Z2))/(2*Z2));
    rmax=inf(((1-Z0-Z1) + sqrt((1-Z0-Z1)^2-4*Y0*Z2))/(2*Z2));
    if rmin<rmax 
      success=1;
      I=[rmin rmax];
      disp(['success: rmin = ',num2str(rmin)])
    else
      disp('failure: rmin > rmax')
    end
  else
    disp('failure: discriminant is negative')  
  end
else
    disp('failure: linear term is positive')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Proof of the eigenvalue problem %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

h = infsup(-rmin,rmin);

tx = [bx(1) + h; bx(2) + h];

A = Dphi_vector_field(tx);

lambda = -0.187256681090725;

v = [-0.581870201791746;-0.813281666009282];

[tlambda,tv] = verifyeig(A,lambda,v);

%save initial_data_ex2 tx tlambda tv
%clear