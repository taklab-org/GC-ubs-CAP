function x = newton_vector_field(x)

tol = 5e-16; %% tolerance for Newton's method

f = phi_vector_field(x);

nf = norm(f,1);
display(['At the beginning ||f|| = ',num2str(nf)])

k=0;

while (k<=100) && (nf > tol)
    Df = Dphi_vector_field(x);
    x = x - Df\f;
    f = phi_vector_field(x);
    nf = norm(f);
    display(['||f|| = ',num2str(nf),', ||Df^(-1)|| = ',num2str(norm(inv(Df),1))])
    k = k+1;                            
end

end