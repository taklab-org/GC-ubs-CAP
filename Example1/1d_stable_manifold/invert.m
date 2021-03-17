function a = invert(b)

% Given a sequence of Taylor coefficients b = (b_n)_{n>=0} 
% representing the function v(t) = \sum_{n\ge0} b_n t^n, 
% this code computes (if possible) the Taylor coefficients a = (a_n)_{n>=0} 
% representing the function u(t)=1/v(t) = \sum_{n\ge0} a_n t^n

tol = 5e-14; %% tolerance for Newton's method

N = length(b)-1 ;

a = [1/b(1);zeros(N,1)];

f = f_invert(a,b);

nf = norm(f,1);
display(['At the beginning ||f|| = ',num2str(nf)])

k=0;
Df = Df_invert(b);

while (k<=100) && (nf > tol)
    a = a - Df\f;
    f = f_invert(a,b);
    nf = norm(f);
    display(['||f|| = ',num2str(nf),', ||Df^(-1)|| = ',num2str(norm(inv(Df),1))])
    k = k+1;                            
end

end

function f = f_invert(a,b)

N = length(a)-1 ;
ab = cauchy([a b]);
delta = [1;zeros(N,1)];
f = mid(ab - delta) ;

end

function Df = Df_invert(b)

N = length(b)-1 ;
Df = zeros(N+1);

for n = 0:N
    j = (0:n);
    k = n-j+1;
    Df(n+1,j+1) = b(k);
end

end


