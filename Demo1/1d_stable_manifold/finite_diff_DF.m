function DF = finite_diff_DF(a,par)

h = 1e-6;
m = length(a);

DF = zeros(m);
E = eye(m);

for j = 1:m
    ah = a + h*E(:,j);
    DF(:,j)=(f_stable_manifold(ah,par) - f_stable_manifold(a,par))/h;
end

end