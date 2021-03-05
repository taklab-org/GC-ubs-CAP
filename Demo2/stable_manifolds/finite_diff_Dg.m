function Dg = finite_diff_Dg(x)

h = 1e-6;

Dg = zeros(3);
E = eye(3);

for j = 1:3
    xh = x + h*E(:,j);
    Dg(:,j)=(g_vector_field(xh) - g_vector_field(x))/h;
end

end