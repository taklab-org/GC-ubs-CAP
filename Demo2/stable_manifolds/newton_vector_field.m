function x = newton_vector_field(x)

tol = 1e-16; %% tolerance for Newton's method

g = g_vector_field(x);

ng = norm(g,inf);
display(['At the beginning ||g|| = ',num2str(ng,inf)])

k=0;

while (k<=100) && (ng > tol)
    Dg = Dg_vector_field(x);
    x = x - Dg\g;
    g = g_vector_field(x);
    ng = norm(g,inf);
    display(['||g|| = ',num2str(ng),', ||Dg^(-1)|| = ',num2str(norm(inv(Dg),inf))])
    k = k+1;                            
end

end