function a = padding(a,n)

N = (length(a)+2)/2; % a = (a1,a2) with ai = (ai_n)_{n=2}^N \in R^{N-1}

if n>0
    a1 = [a(1:N-1);zeros(n,1)]; a2 = [a(N:2*N-2);zeros(n,1)];
else
    a1 = a(1:N-1+n); a2 = a(N:2*N-2+n);
end

a = [a1;a2];

end

