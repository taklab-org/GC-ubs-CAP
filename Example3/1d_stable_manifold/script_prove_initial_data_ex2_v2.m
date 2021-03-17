
clear

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Proof of the fixed point %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

bx = [0.8861081289780320;0.6192579489210105];

ibx = intval(bx);

phi = phi_vector_field(ibx);

Dphi = Dphi_vector_field(bx);

A = intval(inv(Dphi));

Y0 = norm(A*phi,inf);

disp(['Y0 = ',num2str(sup(Y0))])

r_star = 1e-10;

h = infsup(-r_star,r_star);

b1 = bx(1) + h;
b2 = bx(2) + h;

Id = eye(2);

kappa = norm(Id - A*Dphi_vector_field([b1;b2]),inf);

disp(['kappa = ',num2str(sup(kappa))])

Y0 = intval(Y0);

if sup(kappa)>1
    disp('failure: kappa > 1')
    return
else
    rmin = sup(Y0/(1-kappa));
    rmax = r_star;
    disp(['success: rmin = ',num2str(rmin)])
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