function a = newton(a,scaling)

tol = 3e-14; %% tolerance for Newton's method

f = f_1d_stable_manifold_ex3(a,scaling);

nf = norm(f,1);
display(['At the beginning ||f|| = ',num2str(nf)])

k=0;

while (k<=10) && (nf > tol)
    Df = mid(Df_1d_stable_manifold_ex3(a,scaling));
    a = a - Df\f;
    f = f_1d_stable_manifold_ex3(a,scaling);
    nf = norm(f);
    display(['||f|| = ',num2str(nf),', ||Df^(-1)|| = ',num2str(norm(inv(Df),1))])
    k = k+1;                            
end

end