function Df = finite_diff_DF(a,scaling)

h = 1e-6;
m = length(a);

Df = zeros(m);
E = eye(m);

for j = 1:m
    ah = a + h*E(:,j);
    Df(:,j)=(f_1d_stable_manifold_ex3(ah,scaling) - f_1d_stable_manifold_ex3(a,scaling))/h;
end

end