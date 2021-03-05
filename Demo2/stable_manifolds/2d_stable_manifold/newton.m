function a = newton(a,scaling,fp)

tol = 2e-14; %% tolerance for Newton's method

f = f_2d_stable_manifold_ex3(a,scaling,fp);

nf = norm(f,1);
display(['At the beginning ||f|| = ',num2str(nf)])

k=0;

while (k<=10) && (nf > tol)
    Df = mid(Df_2d_stable_manifold_ex3(a,scaling,fp));
    a = a - Df\f;
    f = f_2d_stable_manifold_ex3(a,scaling,fp);
    nf = norm(f);
    display(['||f|| = ',num2str(nf),', ||Df^(-1)|| = ',num2str(norm(inv(Df),1))])
    k = k+1;                            
end

end