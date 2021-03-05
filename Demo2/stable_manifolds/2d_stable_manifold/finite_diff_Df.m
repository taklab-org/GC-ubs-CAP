function Df = finite_diff_Df(a,scaling,fp)

h = 1e-6;

m = length(a);

Df = zeros(m);
E = eye(m);

for j = 1:m
    ah = a + h*E(:,j);
    Df(:,j)=(f_2d_stable_manifold_ex3(ah,scaling,fp) - f_2d_stable_manifold_ex3(a,scaling,fp))/h;
end

end