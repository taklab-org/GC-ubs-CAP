function a = newton(a,par)

tol = 5e-14; %% tolerance for Newton's method

f = f_stable_manifold(a,par);

nf = norm(f,1);
display(['At the beginning ||f|| = ',num2str(nf)])

k=0;

while (k<=100) && (nf > tol)
    Df = mid(Df_stable_manifold(a,par));
    a = a - Df\f;
    f = f_stable_manifold(a,par);
    nf = norm(f);
    display(['||f|| = ',num2str(nf),', ||Df^(-1)|| = ',num2str(norm(inv(Df),1))])
    k = k+1;                            
end

end