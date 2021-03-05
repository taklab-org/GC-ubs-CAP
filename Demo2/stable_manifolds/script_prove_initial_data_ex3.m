%clear

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Proof of the fixed point %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load pts

% Choose the fixed point
point = 1;

if point == 0
    bx = p0;
    v = v0;
    lambda = lbda0;
end

if point == 1
    bx = p1;
    V = V1;
    LBDA = LBDA1;
end

if point == 2
    bx = p2;
    V = V2;
    LBDA = LBDA2;
end

ibx = intval(bx);

g = g_vector_field(ibx);

Dg = Dg_vector_field(bx);

A = intval(inv(Dg));

Y0 = norm(A*g,inf);

disp(['Y0 = ',num2str(sup(Y0))])

r_star = 1e-10;

h = infsup(-r_star,r_star);

b1 = bx(1) + h;
b2 = bx(2) + h;
b3 = bx(3) + h;

Id = eye(3);

kappa = norm(Id - A*Dg_vector_field([b1;b2;b3]),inf);

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

tx = [bx(1) + h; bx(2) + h; bx(3) + h];

A = Dg_vector_field(tx);

if point == 0
    [tlambda0,tv0] = verifyeig(A,lbda0,v0);
    ip0 = tx;
end

if point == 1
    [tlambda1_1,tv1_1] = verifyeig(A,LBDA1(1),V1(:,1));
    [tlambda1_2,tv1_2] = verifyeig(A,LBDA1(2),V1(:,2));
    ip1 = tx;
end

if point == 2
    [tlambda2_1,tv2_1] = verifyeig(A,LBDA2(1),V2(:,1));
    [tlambda2_2,tv2_2] = verifyeig(A,LBDA2(2),V2(:,2));
    ip2 = tx;
end

%save initial_data_ex2 tx tlambda tv
%clear