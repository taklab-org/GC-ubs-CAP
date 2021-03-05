function str = dispintval(b)
b_inf = num2str(inf(b),16);
b_sup = num2str(sup(b),16);
n = diffbit(b_inf,b_sup);
str = [b_inf(1:n-1),'_{',b_inf(n:end),'}^{',b_sup(n:end),'}'];